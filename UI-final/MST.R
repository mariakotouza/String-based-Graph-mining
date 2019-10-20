


MST<-function(x,algorithm="prim"){
forest=components(x,mode="weak")

    no=forest$no
    member=forest$membership
    overall=c()
    for (i in 1:no)
          {ig=induced_subgraph(x,which(member==i),impl="auto")
          plot(ig)
          V(ig)$id2=V(ig)$id
          V(ig)$id=1:gorder(ig)
          vertices=as_data_frame(ig,what="vertices")
          edges=as_data_frame(ig,what="edges")
          
          edges_2=getMinimumSpanningTree(vertices$id,as.matrix(edges[,c("from","to","weight")]),algorithm=algorithm,show.data = FALSE,show.graph=FALSE)$tree.arcs
          edges_2[,"ept1"]=V(ig)$id2[edges_2[,"ept1"]]
          edges_2[,"ept2"]=V(ig)$id2[edges_2[,"ept2"]]
          
          overall=rbind(overall,edges_2)}
      
    overall=as.data.frame(overall)
    colnames(overall)<-c("from","to","weight")
    return(overall)}