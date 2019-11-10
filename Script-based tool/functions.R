# Notes: The weight matrix containes the normalized distances, which means that 0 corresponds to zero distance - we use 0.00001 instead of 0 in this case.
#       During the graph sparsification, we set w[i,j] to zero if there is no connection between i and j.


read_data <- function(datapath, sep="\t", header=T){
  data <<- read.csv(file=as.character(datapath), sep=sep, header = header)#[1:50,]
  #data <- data[sample(1:nrow(data),500, replace = FALSE),]
  #sim=stringdistances(seq=as.character(data[,seqSelect]),algo=simSelect)
  values$indexes <<- 1:nrow(data)
  values$uni <<- NULL
}

calculate_distance <- function(seqSelect, simSelect){
  cat(paste0("calculate_distance","\t"), file=logFile, append=TRUE)
  cat(paste0(seqSelect, ",", simSelect,"\t"), file=logFile, append=TRUE)
  cat(paste0(nrow(data),"\t"), file=logFile, append=TRUE)
  cat(paste0(ncol(data),"\t"), file=logFile, append=TRUE)
  cat(paste0(Sys.time(),"\t"), file=logFile, append=TRUE)
  
  sim=stringdistances(data[,seqSelect],algo=simSelect)
  sim[sim==0]=0.00001
  #shinyalert("Distances Calculated", type = "success")
  values$sim <<- as.matrix(sim)[values$indexes,values$indexes]
  
  write.table(values$sim, paste0(output_folder, "/", "similarity_matrix_", simSelect,".txt"), sep = "\t")
  # log time end and memory used 
  cat(paste0(Sys.time(),"\t"), file=logFile, append=TRUE)
  cat(pryr::mem_used(), file=logFile, append=TRUE, sep = "\n")
}

graph_construction <- function(clusterId=F, seqSelect, graph.type = "exp",alpha,alpha1,alpha2,k,
                               graph_thr, delete_edges_of_diff_classes){
  if (is.null(values$sim )) {return()}
  templist=list(data[values$indexes,],values$sim)
  
  if (clusterId){
    
    if (is.null(values$uni)){ 
      uni=unique(templist[[1]][,c(seqSelect,"cluster_id")])
      ind=1
      uniset=which(templist[[1]][,seqSelect]==uni[1,1] & templist[[1]][,"cluster_id"]==uni[1,2])
      univec=paste(rownames(templist[[1]])[uniset[-1]],collapse=",")
      for (i in 2:nrow(uni))
      { floorindex=ind[i-1]+1
      uniset=which(templist[[1]][floorindex:nrow(templist[[1]]),seqSelect]==uni[i,1] & templist[[1]][floorindex:nrow(templist[[1]]),"cluster_id"]==uni[i,2])
      univec=c(univec,paste(rownames(templist[[1]])[ind[i-1]+uniset[-1]],collapse=","))
      ind=c(ind,ind[i-1]+uniset[1])
      }
      values$uni <<- ind 
      values$uniframe <<- data.frame("Same Combination"=univec,row.names=rownames(templist[[1]])[ind])
    }
    else{ 
      ind=values$uni
    }
    
    templist[[1]]=templist[[1]][ind,]
    templist[[2]]=templist[[2]][ind,ind]
  }
  
  if (delete_edges_of_diff_classes){
    actual_classes <- data$V2
    labels <- unique(actual_classes)[which(!is.na(unique(actual_classes)))]
    for (cl in 1:length(labels)){
      i1 <- which(actual_classes == labels[cl])
      i2 <- which((actual_classes != labels[cl]) & (!is.na(actual_classes)))
      if (length(i2)>0){
        templist[[2]][i1,i2] <- 0.99999
        templist[[2]][i2,i1] <- 0.99999
      }
    }
  }
  d <- templist[[2]]
  print(d[1:10,1:10])
  
  #templist[[2]][templist[[2]]>graph_thr]=0

  # Graph sparsification 
  if ((graph.type == "knn") | (graph.type == "mknn"))
  {
    index <- sapply(1:nrow(d),function(i)
    {
      return(sort(d[,i],index.return = TRUE)$ix[1:k])
    })
    w <- matrix(0,ncol = nrow(d),nrow = nrow(d))
    for (i in 1:nrow(d))
    {
      w[index[,i],i] <- d[index[,i],i]
    }
  }
  if (graph.type == "mknn"){
    for (i in 1:nrow(d)){
      for (j in 1:nrow(d)){
        if (w[i,j]!=w[j,i]){
          w[i,j] <- 0
          w[j,i] <- 0
        }
      }
    }
  }
  
  if (graph.type == "enn")
    w <- ifelse(d < graph_thr,d,0)
  if (graph.type == "tanh")
    w <- (tanh(alpha1 * (d - alpha2)) + 1) / 2
  if (graph.type == "exp")
    w <- exp(-d ^ 2 / alpha ^ 2)
  
  templist[[2]] <- w
  
  ig<- graph.adjacency(templist[[2]], mode="undirected", weighted=TRUE)
  
  ig=delete.edges(ig,which(E(ig)$weight==0))
  
  vertex.attributes(ig)<-as.list(templist[[1]])
  
  V(ig)$value=log10(V(ig)$freq_cluster_id/100)#V(ig)$freq_cluster_id^(1/3.5)
  
  myFun<-function(x)
  {
    if (length(x)==0)
      return (NA)
    else
      return (x[1])
  }
  
  if (domain == "AA"){
    V(ig)$V.GENE=as.character(unlist(lapply(strsplit(as.character(V(ig)$Summary.V.GENE.and.allele),"*",fixed=TRUE), myFun)))
    V(ig)$J.GENE=as.character(unlist(lapply(strsplit(as.character(V(ig)$Summary.J.GENE.and.allele),"*",fixed=TRUE), myFun)))
    V(ig)$D.GENE=as.character(unlist(lapply(strsplit(as.character(V(ig)$Summary.D.GENE.and.allele),"*",fixed=TRUE), myFun)))
  }
  
  V(ig)$color.background=colors[visualiseGenes(data=V(ig)$dataName)$label]
  V(ig)$color.border=V(ig)$color.background
  #!!!!!!!!!!!!!!!!
  V(ig)$title<-paste("Sequence frequence cluster_id: ",V(ig)$freq_cluster_id,"<br>V.Gene: ","<br>CDR3: ")
  
  m=(1-E(ig)$weight-min(1-E(ig)$weight))/(max(1-E(ig)$weight)-min(1-E(ig)$weight))
  
  
  E(ig)$width <- m*2.5
  E(ig)$title<-paste("Edge weight: ",E(ig)$weight)
  
  V(ig)$id=1:gorder(ig)
  V(ig)$label<-rownames(templist[[1]])
  
  forest=components(ig,mode="weak")
  
  x=forest[[1]]
  y=table(x)
  z=y[match(x,names(y))]
  x[z<15]=NA
  uni=cbind(unique(x),y[unique(x)])
  
  
  x=match(x,sort(unique(x),na.last = TRUE))
  values$forest <<- x
  values$ig <<- ig
  values$sim <<- templist[[2]]
  
  
  y=table(values$forest)
  choices=as.numeric(rownames(y))
  names(choices)=paste("Component" ,rownames(y),", ",y," nodes")
  names(choices)[length(choices)]=paste(names(choices)[length(choices)]," (Solo Nodes)")
  #updateSelectInput(session,"componentSelect",choices=choices,selected=1)
  clusterValues$member=NULL
  clusterValues$membership=vector(mode="list",length=7)
  clusterValues$comm1=NULL
  clusterValues$comm2=NULL
  
  names(clusterValues$membership)=c("louvain","fast_greedy","label_propagation","leading_eigenvalue","walktrap","edge_betweenness","hierarchical")
  #updateSelectInput(session,"clusterSelect",selected=" ")
  mstValues$edges=NULL
  
  return(ig)
}

select_component <- function(componentSelect){
  if (is.null(values$ig)) {return()}
  values$ig <<- induced_subgraph(values$ig,which(values$forest==componentSelect),impl="auto")
  V(values$ig)$id <<- 1:gorder(values$ig)
  values$forest <<- c()
  values$forest[1:gorder(values$ig)] <<- 1
  #updateSelectInput(session,"componentSelect",choices=list("1"=1),selected=1)
  clusterValues$member <<- NULL
  clusterValues$membership=vector(mode="list",length=7)
  clusterValues$comm1 <<- NULL
  clusterValues$comm2 <<- NULL
  names(clusterValues$membership) <- c("louvain","fast_greedy","label_propagation","leading_eigenvalue","walktrap","edge_betweenness","hierarchical")
  #updateSelectInput(session,"clusterSelect",selected=" ")
  mstValues$edges <<- NULL
  clusterValues <<- clusterValues
}

adj_table <- function(){
  igtemp=values$ig
  x=adjacent_vertices(igtemp,1:gorder(igtemp))
  
  myFun<-function(x,ig){
    x=as.numeric(x)
    x=V(ig)$label[x]
    neigh=paste(x,collapse=",")
    return (neigh)
  }
  
  adj_table <- data.frame("Neighbours"=unlist(lapply(x,myFun,igtemp)),row.names=V(igtemp)$label)
  write.table(adj_table, paste0(output_folder,"/","adj_table.txt"), row.names = F, sep = "\t")
}

visualize_network <- function(graph_thr) {
  #xx=input$visGraph
  igtemp=isolate(values$ig)
  if (is.null(igtemp)) {return()}
  
  coords=layout_with_stress(igtemp)

  data_vertices=get.data.frame(igtemp,what=c("vertices"))[,c("value","label","color.border","color.background","id")]
  data_edges=as_data_frame(igtemp,what=c("edges"))
  
  graph<-visNetwork(cbind(data_vertices,"title"=paste(V(igtemp)$title,"<br>Component ",isolate(values$forest))), data_edges) %>%
    visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
    visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
    #visPhysics(solver="repulsion") %>%
    visInteraction(multiselect = TRUE)%>%
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE) %>%
    #visIgraphLayout(layout = "layout_with_fr",smooth=FALSE,physics = FALSE) %>%
    visLegend(addNodes = data.frame()) 
  
  #png(filename = paste0("regular_graph_thr_", graph_thr,".png"))
  #print(graph)
  #dev.off()
  
  html_name <- tempfile(fileext = ".html")
  visSave(graph, html_name)
  webshot(html_name, zoom = 2, file = paste0(output_folder,"/","regular_graph_thr_", graph_thr,".png"))
}

filter_data <- function(excludeButton, includeButton, select1, slider1, text1, ReverseButton){
  if(excludeButton[[1]]==0 & includeButton[[1]]==0) {return()}
  tempdata=data
  
  if (isolate(values$flagEX)!=excludeButton[[1]]){
    ie="E"
    values$flagEX=excludeButton[[1]]
  }else {ie="I"}
  
  selected=select1
  if (is.numeric(tempdata[,selected])){ 
    x=data.frame("Columns"=selected,"Keys"=paste0(slider1[1],",",slider1[2]),"I/E"=ie)
  }else{
    x=data.frame("Columns"=selected,"Keys"=text1,"I/E"=ie)
  }
  
  
  values$filterdf <<- rbind(values$filterdf,x)
  
  if (nrow(values$filterdf)==0) {return()}
  tempdata=isolate(data)
  filterg=values$filterdf
  exc=filterg[filterg[,3]=="E",1:2]
  inc=filterg[filterg[,3]=="I",1:2]
  filterinc=list()
  filterexc=list()
  
  if ((nrow(exc))!=0){
    for (i in 1:nrow(exc)){
      filterexc[[i]]=c(as.character(exc[i,1]),unlist(strsplit(as.character(exc[i,2]),",")))
    }}
  
  if ((nrow(inc))!=0){
    for (i in 1:nrow(inc)){
      filterinc[[i]]=c(as.character(inc[i,1]),unlist(strsplit(as.character(inc[i,2]),",")))
    }}  
  
  
  pointers1=includeInGraph(tempdata,filterinc)
  if (length(pointers1)!=0) {
    pointers2=excludeFromGraph(tempdata[pointers1,],filterexc)
    if (length(pointers2)!=0){
      values$indexes <<- pointers1[pointers2]
    }else {
      values$indexes <<- NULL
    }
  }else{
    values$indexes <<- NULL
  }
  
  if (ReverseButton==TRUE){values$indexes <<- setdiff(1:nrow(data),values$indexes)}
}

download_graph <- function(){
  fileName<- paste0(output_folder,"/","edges.txt")
  write.table(as_data_frame(values$ig)[,c("from","to","weight")], fileName, row.names = FALSE, sep = "\t")
  files=fileName
  fileName<-paste0(output_folder,"/","vertices.txt")
  write.table(V(values$ig)$label, fileName, row.names = FALSE, sep = "\t")
}

neigh_table <- function(){
  neighbors
}

create_mst <- function(mstAlgoSelect) {
  if (is.null(values$ig)) {return()}
  igtemp=values$ig
  if (mstAlgoSelect=="Kruskal"){
    edges=MST(igtemp,mstAlgoSelect)
  }else{
    edges=as_data_frame(mst(igtemp,algorithm = "prim"))
  }
  #a=mstClustering(values$ig)
  ig=graph_from_data_frame(edges,directed=FALSE,vertices=as_data_frame(igtemp,what="vertices")[,c("id")])
  
  degree=centr_degree(ig)
  centroids=which(degree$res>1/15*gorder(ig))
  
  
  n=c()
  n[gorder(ig)+1]=1
  n[centroids]="triangle"
  mstValues$shape <<- n[1:gorder(ig)]

  #V(mstig)$value=V(values$ig)$freq_cluster_id^(1/3.5)*1
  
  m=(1-edges$weight-min(1-edges$weight))/(max(1-edges$weight)-min(1-edges$weight))
  edges$width <- m*2.5
  
  #V(mstig)$id=1:gorder(mstig)
  mstValues$edges <<- edges
  mstValues$centroids <<- V(igtemp)$label[centroids] #a$centroids]
  mstValues$flag2 <<- NULL
}

adj2_table <- function(){
  req(mstValues$edges,isolate(values$ig))
  edges=mstValues$edges
  igtemp=isolate(values$ig)
  ig=graph_from_data_frame(edges,directed=FALSE,vertices=as_data_frame(igtemp,what="vertices")[,c("id","label")])
  x=adjacent_vertices(ig,1:gorder(ig))
  
  
  
  
  myFun<-function(x,ig){
    x=as.numeric(x)
    x=V(ig)$label[x]
    neigh=paste(x,collapse=",")
    return (neigh)
  }
  
  adj2_table <- data.frame("Neighbours"=unlist(lapply(x,myFun,ig)),row.names=V(ig)$label)
  write.table(adj2_table, paste0(output_folder,"/","adj_table_mst.txt"), row.names = F, sep = "\t")
}

mst_construction <- function(colormst, bordermst, clusterMST){
  
  #xx=input$visMST
  #isolate({
  igtemp=values$ig
  edges=mstValues$edges
  if (is.null(edges) | is.null(igtemp)) {return()}
  
  ig=graph_from_data_frame(edges,directed=FALSE,vertices=as_data_frame(igtemp,what="vertices")[,c("id","label","value","title")])
  
  degree=centr_degree(ig)
  centroids=which(degree$res>1/15*gorder(ig))
  coords=layout_with_stress(ig)
  #coords=layout_with_kk(igtemp,coords = coords)
  
  V(ig)$id=V(ig)$name
  delete_vertex_attr(ig,"name")
  
  
  if ((clusterMST) & !is.null(clusterValues$member) ){   
    cbg=colors[clusterValues$member]
    mstValues$legenddf=data.frame("label"=unique(clusterValues$member),"color"=(colors[unique(clusterValues$member)]))
  }else{  
    if (colormst=="Default"){
      cbg="black"
      mstValues$legenddf=NULL
    }else{
      temp=visualiseGenes(get.vertex.attribute(igtemp,colormst))
      cbg=colors[temp$label]
      mstValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])}
  }#unique(temp$label)
  
  if (bordermst=="Default"){
    cbor="black" #colors[a$clusters]
    mstValues$legenddf2=NULL
  }else{
    temp=visualiseGenes(get.vertex.attribute(igtemp,bordermst))
    cbor=colors[temp$label]
    mstValues$legenddf2=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])
  }
  
  mstValues <<- mstValues
  
  write.table(mstValues$edges[,c("from","to")], paste0(output_folder,"/","mst.txt"), row.names = FALSE,sep="\t")
  
  result <- list("ig"=ig, "cbg"=cbg, "cbor"=cbor, "coords"=coords)
  return(result)
}

mst_visualize <- function(colormst, bordermst, clusterMST, ig, cbg, cbor, coords){
  
  graph_mst <- visNetwork(cbind(as_data_frame(ig,what="vertices"),"color.background"=cbg,"color.border"=cbor,"shape"=mstValues$shape),cbind(as_data_frame(ig),"title"=paste("Edge weight: ",E(ig)$weight))) %>%
      visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
      visNodes(borderWidth = 2.5) %>%
      visOptions (nodesIdSelection = TRUE,highlightNearest = TRUE)   %>%
      ##visPhysics(solver="forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant=-2000,centralGravity=0.02,springLength=40,springConstant=0.4,damping=1,avoidOverlap=1)) %>%
      visInteraction(multiselect = TRUE)%>%
      visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE) 
    
    html_name <- tempfile(fileext = ".html")
    visSave(graph_mst, html_name)
    webshot(html_name, zoom = 2, file = paste0(output_folder,"/","graph_mst",".png"))
    
    
  #})
  
  if (!is.null(mstValues$legenddf)){
    mstValues$legenddf=mstValues$legenddf[order(mstValues$legenddf$label),]
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("bottomleft", legend =mstValues$legenddf$label, pch=16, pt.cex=3, cex=1.5, bty='n',
         col = as.character(mstValues$legenddf$color))
    mtext("Background",side=2, at=0.2, cex=2)
  }
  
  if (!is.null(mstValues$legenddf2)){
    mstValues$legenddf2=mstValues$legenddf2[order(mstValues$legenddf2$label),]
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("bottomleft", legend =mstValues$legenddf2$label, pch=16, pt.cex=3, cex=1.5, bty='n',
         col = as.character(mstValues$legenddf2$color))
    mtext("Border",side=2, at=0.2, cex=2)
  }
 
}

mst_with_clustering_colors <- function(colormst, bordermst, clusterColor, centralColor, clusterSelect2, clusterMST){
  x=colormst
  y=bordermst
  z=clusterMST
  
  if ( is.null(isolate(mstValues$flag2)) | is.null(isolate({values$ig}))){return()}
  tempdata=as_data_frame(isolate({values$ig}),what="vertices")
 
  
  if (!is.null(clusterValues$member) & clusterMST){
    bgcolor=colors[clusterValues$member]
    mstValues$legenddf=data.frame("label"=unique(clusterValues$member),"color"=colors[unique(clusterValues$member)])
  }else{  
    if (x=="Default"){
      bgcolor="black"#colors[mstValues$clusters]
      mstValues$legenddf=NULL
      #data.frame("label"=unique(mstValues$clusters),"color"=unique(colors[mstValues$clusters]))
    }else{
      temp=visualiseGenes(tempdata[,colormst])
      bgcolor=colors[temp$label]
      mstValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])
    }
  }
  
  if (y=="Default"){
    bordcolor="black" #colors[mstValues$clusters]
    mstValues$legenddf2=NULL 
    #data.frame("label"=unique(mstValues$clusters),"color"=unique(colors[mstValues$clusters]))
  }else{
    temp=visualiseGenes(tempdata[,bordermst])
    bordcolor=colors[temp$label]
    mstValues$legenddf2=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])
  }
  
  
  visNetwork(data.frame(id=1:nrow(tempdata),"color.background"=bgcolor,"color.border"=bordcolor,"shape"=mstValues$shape),cbind(as_data_frame(ig),"title"=paste("Edge weight: ",E(ig)$weight))) %>%
    visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
    visNodes(borderWidth = 2.5) %>%
    visOptions (nodesIdSelection = TRUE,highlightNearest = TRUE)   %>%
    ##visPhysics(solver="forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant=-2000,centralGravity=0.02,springLength=40,springConstant=0.4,damping=1,avoidOverlap=1)) %>%
    visInteraction(multiselect = TRUE)%>%
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE) 
  
}

calculate_centrality <- function(fastCButton=F, centralSelect, clusterCentral, centralColor){
  if (is.null(values$ig)){return()}
  
  cat(paste0("calculate_centralities","\t"), file=logFile, append=TRUE)
  cat(paste0(fastCButton, "=", fastCButton, ",", centralSelect,"\t"), file=logFile, append=TRUE)
  cat(paste0("\t"), file=logFile, append=TRUE)
  cat(paste0("\t"), file=logFile, append=TRUE)
  cat(paste0(Sys.time(),"\t"), file=logFile, append=TRUE)
  
  igtemp=giant_component_extract(values$ig)[[1]]
  
  prop=proper_centralities(igtemp)
  print("centralities")
  if (fastCButton){
    print(system.time({
      central=calculate_centralities(igtemp,include=prop[c(3,10,26)])
      E(igtemp)$weight=1-E(igtemp)$weight
      central=c(central,calculate_centralities(igtemp,include=prop[15]))
    }))
    b=c("Average.Distance")
  }else { 
    print(system.time({
      central=calculate_centralities(igtemp,include=prop[c(-2,-5,-6,-8,-15,-19:-20,-21,-24,-27,-31:-33,-38:-39,-43,-50)])
      E(igtemp)$weight=1-E(igtemp)$weight
      central=c(central,calculate_centralities(igtemp,include=prop[c(15,19,20,24,50)]))
    }))
    b=c("Average.Distance","Local.Bridging.Centrality","Wiener.Index.Centrality")
  }
  
  central=central[lengths(central)!=0]
  x=as.data.frame(central,col.names = names(central),row.names=V(igtemp)$label)
  
  print("ranks")  
  print(system.time({
    central_ranked=as.data.frame(apply(x,2,sort))
    colnames(central_ranked)=colnames(x)
    
    
    #revert the order of indices 
    central_ranked[,-match(b,colnames(central_ranked))]=central_ranked[nrow(central_ranked):1,-match(b,colnames(central_ranked))]
    
    ranks<-function(i,central_ranked,x){
      positions_first=match(x[,i],central_ranked[,i])
      positions_last=match(x[,i],central_ranked[nrow(x):1,i])
      return (positions_first-positions_last+nrow(x)+1)/2}
    
  }))
  
  positions<-as.data.frame(lapply(1:length(x),ranks,central_ranked=central_ranked,x=x))
  
  colnames(positions)=colnames(x)
  rownames(positions)=rownames(x)
  
  centralValues$centralities=x
  centralValues$rankings=positions
  centralValues$summary=as.data.frame(t(apply(positions,1,summary)),rownames=rownames(x))
  
  #Network 
  igtemp=giant_component_extract(isolate(values$ig))[[1]]
  if (is.null(igtemp) | is.null(centralValues$centralities)) {return()}
  
  
  bc <- centralValues$centralities[,centralSelect]
  if (centralSelect=="Average.Distance") {bc=1/bc}
  
  coords=layout_with_centrality(igtemp,cent=bc,iter=500,tol=1e-04)
  # p1 <- ggraph(igtemp,layout = "centrality", cent = bc,iter = 500, tol = 1e-04) +
  #   draw_circle(use = "cent") +
  #   annotate_circle(bc,format="",pos="bottom") +
  #   geom_edge_link(edge_color="black",edge_width=0.3)+
  #   geom_node_point(aes(fill=as.factor(Faction)),size=2,shape=21)+
  #   scale_fill_manual(values=c("#8B2323", "#EEAD0E"))+
  #   theme_graph()+
  #   theme(legend.position = "none")+
  #   coord_fixed()+
  #   labs(title="betweenness centrality")
  # 
  # coords=cbind(x=p1$data$x,y=p1$data$y)
  
  
  
  i=c()
  i[gorder(igtemp)+1]=""
  i[coords[,1]==0 & coords[,2]==0]="triangle"
  
  
  
  temp=visualiseGenes(get.vertex.attribute(igtemp,"dataName"))
  centralValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])
  #updateSelectInput(session,"centralColor",selected="dataName")
  #updateCheckboxInput(session,"clusterCentral",value=FALSE)
  
  data_vertices=cbind(as_data_frame(igtemp,what=c("vertices"))[,c("value","label","title")],shape=i[1:gorder(igtemp)],color=colors[temp$label])
  data_vertices$id=rownames(data_vertices)
  data_edges=data.frame(from=1,to=2,weight=1)
  
  visNetwork(data_vertices, data_edges) %>%
    visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
    visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
    #visPhysics(solver="repulsion") %>%
    visInteraction(multiselect = TRUE)%>%
    visEdges(color="white")%>%
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE)
  
  #Colors
  x=clusterCentral
  y=centralColor
  
  if (is.null(isolate(values$ig))){return()}
  igtemp=giant_component_extract(isolate(values$ig))[[1]]
  
  if (!is.null(clusterValues$member) && clusterCentral){
    tempcolor=colors[clusterValues$member[get.vertex.attribute(igtemp,"id")]]
    centralValues$legenddf=NULL
  }else{
    temp=visualiseGenes(get.vertex.attribute(igtemp,centralColor))
    tempcolor=colors[temp$label]
    centralValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])
  }
  
  #visNetworkProxy("centralnetwork") %>%
    #visUpdateNodes(data.frame(id=1:gorder(igtemp),color=tempcolor))
  
  visNetwork(data_vertices, data_edges) %>%
    visEdges(color=list(color=tempcolor,highlight="red"),labelHighlightBold=FALSE)%>%
    visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
    #visPhysics(solver="repulsion") %>%
    visInteraction(multiselect = TRUE)%>%
    visEdges(color="white")%>%
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE)
  
  if (is.null(centralValues$legenddf)){return()}
  centralValues$legenddf=centralValues$legenddf[order(centralValues$legenddf$label),]
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("bottomleft", legend =centralValues$legenddf$label, pch=16, pt.cex=3, cex=1.5, bty='n',
         col = as.character(centralValues$legenddf$color))
  mtext("Nodes Color",side=2, at=0.2, cex=2)
  
  centralValues <<- centralValues
  
  #Download
  fileName<- paste0(output_folder,"/","centralities.txt")
  write.table(cbind(id=row.names(centralValues$centralities), centralValues$centralities), fileName, row.names = FALSE,sep="\t")
  
  fileName<-paste0(output_folder,"/","rankings.txt")
  write.table(cbind(id=row.names(centralValues$rankings),centralValues$rankings), fileName, row.names = FALSE,sep="\t")

  fileName<-paste0(output_folder,"/","summary.txt")
  write.table(cbind(id=row.names(centralValues$summary),centralValues$summary), fileName, row.names = FALSE,sep="\t")
  
  # log time end and memory used 
  cat(paste0(Sys.time(),"\t"), file=logFile, append=TRUE)
  cat(pryr::mem_used(), file=logFile, append=TRUE, sep = "\n")
}

clustering <- function(clusterSelect){
  success <- F
  tryCatch({
  igtemp=values$ig
  if (is.null(igtemp) ){return()}
  
  cat(paste0("clustering","\t"), file=logFile, append=TRUE)
  cat(paste0(clusterSelect, "=", clusterSelect,"\t"), file=logFile, append=TRUE)
  cat(paste0("\t"), file=logFile, append=TRUE)
  cat(paste0("\t"), file=logFile, append=TRUE)
  cat(paste0(Sys.time(),"\t"), file=logFile, append=TRUE)
  
  t1 <- Sys.time()
  if (clusterSelect==" ") {
    clusterValues$member=NULL
  }else{    
    if (is.null(clusterValues$membership[[clusterSelect]])){
      #id <<- showNotification(paste("Clustering Calculation..."), duration = NULL)
      print("cluster")
      print(system.time({model=communities_general(igtemp,algorithm = clusterSelect)}))
      clusterValues$membership[[clusterSelect]]=model$membership
      if (!is.null(id)){
        removeNotification(id)
        id <<- NULL
      }
    }
    clusterValues$member=clusterValues$membership[[clusterSelect]]
  }
  
  clusterValues <<- clusterValues
  
  file <- paste0(output_folder,"/","clustering_", clusterSelect, ".txt")
  write.table(clusterValues$member,file, row.names = FALSE,sep="\t")
  
  # log time end and memory used 
  cat(paste0(Sys.time(),"\t"), file=logFile, append=TRUE)
  cat(pryr::mem_used(), file=logFile, append=TRUE, sep = "\n")
  t2 <- Sys.time()
  print(as.numeric(t2-t1))
  success <- T
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  if (success){
    return(as.numeric(t2-t1))
  }else{
    return(-1)
  }
}

visualize_clustered_graph <- function(clusterColor, clusterSelect){
  #xx=input$visCluster
  igtemp=isolate(values$ig)
  membertemp=isolate(clusterValues$membership[[clusterSelect]])
  if (is.null(membertemp) | is.null(igtemp)) {return()}
  #updateSelectInput(session,"clusterColor",selected="Default")
  
  V(igtemp)$grp=membertemp
  print("layout")
  print(system.time({bb <- layout_as_backbone(igtemp,keep=0.4)}))
  V(igtemp)$color=colors[membertemp]
  
  data_vertices=cbind(as_data_frame(igtemp,what=c("vertices")))[,c("value","label","color","title")]
  data_vertices$id=rownames(data_vertices)
  data_edges=as_data_frame(igtemp,what=c("edges"))
  
  
  
  visNetwork(data_vertices, data_edges) %>%
    visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
    visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
    #visPhysics(solver="repulsion") %>%
    visInteraction(multiselect = TRUE)%>%
    visIgraphLayout(layout = "layout.norm", layoutMatrix = cbind(x=bb$xy[,1],y=bb$xy[,2]),smooth=FALSE,physics = FALSE) 
  
  #colors
  if (is.null(values$ig)){return()}
  igtemp=values$ig
  if (clusterColor!="Default"){
    temp=visualiseGenes(get.vertex.attribute(igtemp,clusterColor))
    tempcolor=colors[temp$label]
    clusterValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])
  }else {
    if (!is.null(clusterValues$membership[[clusterSelect]])){
      tempcolor=colors[clusterValues$membership[[clusterSelect]]]
    }else{
      tempcolor=NULL
    }
    clusterValues$legenddf=NULL
  }
  if (is.null(tempcolor)){return()}
  #visNetworkProxy("clusterNetwork") %>%
    #visUpdateNodes(data.frame(id=1:gorder(igtemp),color=tempcolor))
  
  data_vertices$color <- tempcolor
  graph_clustering <- visNetwork(data_vertices, data_edges) %>%
    visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
    visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
    #visPhysics(solver="repulsion") %>%
    visInteraction(multiselect = TRUE)%>%
    visIgraphLayout(layout = "layout.norm", layoutMatrix = cbind(x=bb$xy[,1],y=bb$xy[,2]),smooth=FALSE,physics = FALSE) 
  
  html_name <- tempfile(fileext = ".html")
  visSave(graph_clustering, html_name)
  webshot(html_name, zoom = 2, file = paste0(output_folder,"/","graph_clustering_", clusterSelect,".png"))
  
  #legend
  if (is.null(clusterValues$legenddf)){return()}
  clusterValues$legenddf=clusterValues$legenddf[order(clusterValues$legenddf$label),]
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("bottomleft", legend =clusterValues$legenddf$label, pch=16, pt.cex=3, cex=1.5, bty='n',
         col = as.character(clusterValues$legenddf$color))
  mtext("Nodes Color",side=2, at=0.2, cex=2)
  
  #Clustering metrics
  if (is.null(clusterValues$membership[[clusterSelect]])) {return()}
  igtemp=isolate(values$ig)
  E(igtemp)$weight=1-E(igtemp)$weight
  print("mod")
  print(system.time({
    mod=modularity(igtemp,clusterValues$membership[[clusterSelect]])}))
  
  
  print("con")
  print(system.time({
    con=conductance(igtemp,clusterValues$membership[[clusterSelect]])$conductance}))
  
  print("cov")
  print(system.time({
    cov=coverage(igtemp,clusterValues$membership[[clusterSelect]])}))
  paste0("Modularity:",mod,"  Conductance:",con,"  Coverage:",cov)
}

clustering_heatmap <- function(plotHeat){
  if (is.null(clusterValues$member) | !plotHeat) {return()}
  dis=as.matrix(isolate(values$ig[]))
  maxt=max(apply(dis,1,max))
  dis[dis==0]=maxt+0.1
  diag(dis)=0
  #dis=distances(isolate(values$ig))
  temp=order(clusterValues$member)
  dis=dis[temp,temp]
  print("heatmap")
  print(system.time({
    heatmap.2(dis,col=colorpanel(256,"red","orange","yellow"),scale="none",trace = "none", density.info = "none",dendrogram="none",Rowv = NA, Colv = NA)}))
  #heatmap(dis, Rowv = NA, Colv = NA, col = heat.colors(256), revC = TRUE)
}

silhouette <- function(members, clust_method){
  #if (is.null(clusterValues$member)) {return()}
  dis=as.matrix(isolate(values$ig[]))
  maxt=max(apply(dis,1,max))
  dis[dis==0]=maxt+0.1
  diag(dis)=0
  #dis=distances(isolate(values$ig))
  
  #print(system.time({
    #png(filename = "silhouette_", clust_method, ".png")
    #plot(cluster::silhouette(dist=dis,x=members))
    #dev.off()
  #}))
  
  if ((length(unique(members)) > 1) & (length(unique(members)) < length(members))){
    sil <- cluster::silhouette(dist=dis,x=members)
    ss <- mean(sil[, 3])
    #print(sil_grouped)
  }else{
    ss <- NA
  }
  
  return(ss)
}

clustering_common <- function(clusterSelect1, clusterSelect2){
  igtemp=values$ig
  if (is.null(igtemp)){return()}
  
  clusterValues$comm1 <<- clusterValues$membership[[clusterSelect1]]#communities_general(igtemp,algorithm =input$clusterSelect1)$membership
  
  if (clusterSelect2 %in% c("louvain","fast_greedy","label_propagation","leading_eigenvalue","walktrap","edge_betweenness")){
    clusterValues$comm2 <<- clusterValues$membership[[clusterSelect2]]
    #communities_general(igtemp,algorithm =clusterSelect2)$membership}
  }else{
    clusterValues$comm2 <<- get.vertex.attribute(igtemp,clusterSelect2)
  }
  
  clusterValues <<- clusterValues
  
}

clustering_confusionMatrix1 <- function(clusterSelect1, clusterSelect2){
  validate(
    need(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)), "One of the clusterings hasn't been calculated!")
  )
  
  mat=as.data.frame(t(confusion(clusterValues$comm1,clusterValues$comm2)$mat1))
  mat[,c(1,2)]<-mat[,c(2,1)]
  isolate({colnames(mat)<-c(clusterSelect1,clusterSelect2,"Freq")})
  
  fileName<- paste0(output_folder,"/","clustering_confusionMatrix1.txt")
  write.table(mat, fileName, row.names = FALSE,sep="\t")
  
  mat
}

clustering_confusionMatrix2 <- function(clusterSelect1, clusterSelect2){
  req(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)))
  mat=as.data.frame(t(confusion(clusterValues$comm2,clusterValues$comm1)$mat1))
  mat[,c(1,2)]<-mat[,c(2,1)]
  isolate({colnames(mat)<-c(clusterSelect2,clusterSelect1,"Freq")})
  
  fileName<- paste0(output_folder,"/","clustering_confusionMatrix2.txt")
  write.table(mat, fileName, row.names = FALSE,sep="\t")
  
  mat
}

clustering_metrics <- function() {
  # if (is.null(clusterValues$comm1) | is.null(clusterValues$comm2)) {return()}
  req(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)))
  nmi=NMI(clusterValues$comm1,clusterValues$comm2)
  ari=adj.rand.index(clusterValues$comm1,clusterValues$comm2)
  paste0("NMI:",nmi," ARI:",ari)
}

clustering_interheat <- function(heatSelect=0){
  req(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)))
  
  if (heatSelect==0){ 
    x=NULL
  }else{
    x=as.numeric(heatSelect)
  }
  
  heatmap.2(t(prop.table(table(clusterValues$comm1,clusterValues$comm2),x)), scale="none",trace = "none",
            xlab="Clustering 1",ylab="Clustering 2",main="Probability Heatmap",density.info = "none",col=colorpanel(256,"blue","purple","red"),dendrogram = "none",Rowv = NA, Colv = NA)
  
}

# d: the similarity matrix
# y: thw knows labels
# known.label: the ids of the points that have a known label
sslLabelProp_with_algo_select <- function(d,y,known.label,graph.type = "exp",dist.type = "Euclidean",alpha,alpha1,alpha2,k,epsilon,iter = 1000) {
    all.Obs = nrow(w)
    all.label = unique(y)
    C <- length(unique(y))
    num.known = length(known.label)
    if (num.known != length(y))
      stop("the number of known.label  doesn't accord with that of y")
    num.class = dim(y)[2]
    #d <- proxy::dist(x,x,method = dist.type)
    #d <- matrix(d,ncol = all.Obs)
    if (graph.type == "knn")
    {
      index <- sapply(1:all.Obs,function(i)
      {
        return(sort(d[,i],index.return = TRUE)$ix[1:k])
      })
      w <- matrix(0,ncol = all.Obs,nrow = all.Obs)
      for (i in 1:all.Obs)
      {
        w[index[,i],i] <- d[index[,i],i]
      }
    }
    if (graph.type == "enn")
      w <- ifelse(d < epsilon,d,0)
    if (graph.type == "tanh")
      w <- (tanh(alpha1 * (d - alpha2)) + 1) / 2
    if (graph.type == "exp")
      w <- exp(-d ^ 2 / alpha ^ 2)
    p <- NetPreProc::Prob.norm(w)
    p <- t(p)
    rm(w,d)
    ff <- rep(0,all.Obs)
    ff[known.label] <- y
    for (i in 1:iter)
    {
      ff <- ff %*% p
      ff[known.label] <- y
    }
    return (as.vector(ff))
  }

sslLabelProp <- function(w,y,known.label,iter = 1000,centralities) {
  t1 <- Sys.time()
  all.Obs = nrow(w)
  all.label = unique(y)
  C <- length(unique(y))
  num.known = length(known.label)
  if (num.known != length(y))
    stop("the number of known.label  doesn't accord with that of y")
  num.class = dim(y)[2]
  #d <- proxy::dist(x,x,method = dist.type)
  #d <- matrix(d,ncol = all.Obs)
  
  if (length(centralities) > 1){
    imp_node_ids <- which(centralities > median(centralities))
    for (i in imp_node_ids){
      # for rows
      j <- which(w[i,]!=0)
      if (i==1) print( w[i,j] )
      w[i,j] <- ifelse((w[i,j] - centralities[i])>0, (w[i,j] - centralities[i]), 0.001)
      if (i==1) print( w[i,j] )
      # for cols
      j <- which(w[,i]!=0)
      w[j,i] <- ifelse((w[j,i] - centralities[i])>0, (w[j,i] - centralities[i]), 0.001)
    }
  }
  
  p <- NetPreProc::Prob.norm(w) # The transition probability matrix of network - w contains the distances, not the weights
  p <- t(p) # Transposes the rows and columns of matrices - Now p is the probability matrix
  rm(w)
  ff <- rep(0,all.Obs)
  ff[known.label] <- y

  for (i in 1:iter)
  {
    ff <- ff %*% p # %*% is matrix multiplication. For matrix multiplication, you need an m x n matrix times an n x p matrix
    ff[known.label] <- y
  }
  
  t2 <- Sys.time()
  return (list(f=as.vector(ff), t=(as.numeric(t2-t1))))
}
