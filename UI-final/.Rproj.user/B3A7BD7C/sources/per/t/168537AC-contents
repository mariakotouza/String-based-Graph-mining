library(shiny) 
library(stringdist)
library(igraph)
library(visNetwork)
library(CINNA)
library(DT)

source('C:/Users/Paul Kallin/Desktop/Centralities/stringdistances.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/Centralities/visualisation.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/Centralities/filtered.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/Centralities/mstClustering.R', echo=TRUE)


# use the below options code if you wish to increase the file input limit, in this example file input limit is increased from 5MB to 9MB
# options(shiny.maxRequestSize = 9*1024^2)

shinyServer(function(input,output,session){
  
  values <- reactiveValues(flagEX=0,filterdf=data.frame(),indexes=c())

  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
    list1 <- reactive({
        file1 <- input$file
        if(is.null(file1)){return()} 
        data=read.csv(file=file1$datapath, sep=input$sep, header = input$header)
        data=data[sample(max(dim(data)),50, replace = FALSE),]
        sim=stringdistances(seq=as.character(data[,"IMGT.gapped.AA.sequences.V.D.J.REGION"]))
        list(data,sim)
       })
  
    
    ig<- eventReactive(input$graphButton, {
          templist=list1()
          ig<- graph.adjacency(templist[[2]], mode="undirected", weighted=TRUE)
          ig=delete.edges(ig,which(E(ig)$weight>input$slider))
          vertex.attributes(ig)<-as.list(templist[[1]])
          
          V(ig)$value=V(ig)$freq_cluster_id^(1/3.5)*1
          #V(ig)$size=V(ig)$freq_cluster_id^(1/3.5)*5
          
          #V(ig)$color=visualiseGenes(data=V(ig)$dataName)$label
          
          V(ig)$color.background=colors[visualiseGenes(data=V(ig)$dataName)$label]
          V(ig)$color.border=V(ig)$color.background
          #shapes <- c("circle", "square", "csquare", "rectangle")
          #V(ig)$shape <- shapes[V(ig)$color]
          
          m=(1-E(ig)$weight-min(1-E(ig)$weight))/(max(1-E(ig)$weight)-min(1-E(ig)$weight))
          
          
          E(ig)$width <- m*2.5
          
          V(ig)$id=1:gorder(ig)
          V(ig)$label<-rownames(templist[[1]])
          ig
          })
  
        
        output$dataset<-DT::renderDataTable({
          temp=list1()[[1]]
          if (is.null(temp)){return()}
          if (input$idInput!=0)
             {temp=temp[input$idInput==rownames(temp),]}
         
          temp
          
          })
    
    
    
        output$network <- renderVisNetwork({
          igtemp=ig()
          data_vertices=as_data_frame(igtemp,what=c("vertices"))
          data_edges=as_data_frame(igtemp,what=c("edges"))
          visNetwork(data_vertices, data_edges) %>%
            visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
            visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
            #visPhysics(solver="repulsion") %>%
            visInteraction(multiselect = TRUE)%>%
            visIgraphLayout(layout = "layout_with_fr",smooth=FALSE,physics = FALSE) #%>%
          #visClusteringByColor( colors) 
        })
          
          
          
       output$selectbox1<-renderUI({
         tempdata=list1()
         
         choice=list()
         for (l in colnames(tempdata[[1]]))
            choice=c(choice,l)
         names(choice)=colnames(tempdata[[1]])
         selectInput("select1", label = h3("Select box"), choices = choice)
        })
       
       
       output$selectbox2<-renderUI({
         tempdata=list1()
         
         choice=list()
         for (l in colnames(tempdata[[1]]))
           choice=c(choice,l)
         names(choice)=colnames(tempdata[[1]])
         selectInput("select2", label = h3("Select box"), choices = choice)
       })
       
         
      output$textbox1<-renderUI({
        tempdata=list1()
        selected=input$select1
        if (is.null(selected)) {return()}
        if (is.numeric(tempdata[[1]][,selected]))
           {sliderInput("slider1", label = h3("Slider Range"), min = min(tempdata[[1]][,selected]), max = max(tempdata[[1]][,selected]), value = c(min(tempdata[[1]][,selected]),max(tempdata[[1]][,selected]) ))}
        else
            {textInput("text1","Insert Text")}
        
        
        
        
        
      }) 
      output$textbox2<-renderUI({
        tempdata=list1()
        selected=input$select2
        if (is.null(selected)) {return()}
        if (is.numeric(tempdata[[1]][,selected]))
        {sliderInput("slider2", label = h3("Slider Range"), min = min(tempdata[[1]][,selected]), max = max(tempdata[[1]][,selected]), value = c(min(tempdata[[1]][,selected]),max(tempdata[[1]][,selected]) ))}
        else
        {textInput("text2","Insert Text")}
      })
        
 
      
      

      
      xxchange <- reactive({
        paste(input$excludeButton ,input$includeButton)
      })
     
      
      
      observeEvent(xxchange(),
                              {
                                if(input$excludeButton[[1]]==0 & input$includeButton[[1]]==0) {return()}
                                tempdata=list1()  
                                if (isolate(values$flagEX)!=input$excludeButton[[1]]){
                                    selected=input$select2
                                    
                                    if (is.numeric(tempdata[[1]][,selected])){
                                        x=data.frame("Columns"=selected,"Keys"=paste0(input$slider2[1],",",input$slider2[2]),"I/E"="E")}
                                    else{   
                                        x=data.frame("Columns"=selected,"Keys"=input$text2,"I/E"="E")}
                                    
                                    filterg<<-rbind(filterg,x)
                                    values$flagEX=input$excludeButton[[1]]
                                    
                                }
                                else 
                                  {
                                    selected=input$select1
                                    if (is.numeric(tempdata[[1]][,selected])){ 
                                        x=data.frame("Columns"=selected,"Keys"=paste0(input$slider1[1],",",input$slider1[2]),"I/E"="I")}
                                    else{
                                        x=data.frame("Columns"=selected,"Keys"=input$text1,"I/E"="I")}
                                    
                                    filterg<<-rbind(filterg,x)
                                    }
                                 values$filterdf<-filterg
                                
                              })
      
      
      observeEvent(input$ResetButton,{
                    filterg<<-data.frame(Columns=character(),Keys=character(),"I/E"=character())
                    values$filterdf<-filterg
                    })
      
      
      
      observeEvent(input$FilterButton,{
                  tempdata=isolate(list1()[[1]])
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
                  pointers2=excludeFromGraph(tempdata[pointers1,],filterexc)
                  values$indexes=pointers1[pointers2]
                  
                 },ignoreInit = TRUE)                 
      
      
      
      output$Indexes<-renderText({
        if (is.null(values$indexes)) {return()}
        values$indexes
      })
      
      
      
      output$filtertable<-renderTable({
                          
                          if (is.null(values$filterdf)){return()}
                          values$filterdf 
                              })
      
      
      
      mstValues<-reactiveValues(ig=NULL,centroids=NULL,keyvertices=NULL,clusters=NULL,flag2=NULL)
      
      
      observeEvent(input$mstButton,{
          
          if (is.null(ig())) {return()}
          mstig=mst(ig(),algorithm='prim')
          a=mstClustering(ig())
          
         
          
          if (input$colormst=="Default")
             {V(mstig)$color.background=colors[a$clusters]}
          else
             {V(mstig)$color.background=colors[visualiseGenes(get.vertex.attribute(mstig,input$colormst))$label]}
          
          if (input$bordermst=="Default")
            {V(mstig)$color.border="black"}
          else
            {V(mstig)$color.border=colors[visualiseGenes(get.vertex.attribute(mstig,input$bordermst))$label]}
          
          V(mstig)$value=V(ig())$freq_cluster_id^(1/3.5)*1
          
          m=(1- E(mstig)$weight-min(1-E(mstig)$weight))/(max(1-E(mstig)$weight)-min(1-E(mstig)$weight))
          E(mstig)$width <- m*2.5
          
          V(mstig)$id=1:gorder(mstig)
          mstValues$ig=mstig
          mstValues$centroids=V(mstig)$label[a$centroids]
          mstValues$keyvertices=V(mstig)$label[unique(unlist(a$keyvertices))]
          mstValues$clusters=a$clusters
          mstValues$flag2=NULL
        
      },ignoreInit = TRUE)
      
      
      output$mstnetwork<-renderVisNetwork({
        
        if (is.null(mstValues$ig)) {return()}
        data_vertices=as_data_frame(mstValues$ig,what=c("vertices"))
        data_edges=as_data_frame(mstValues$ig,what=c("edges"))
        mstValues$flag2=c(1)
        visNetwork(data_vertices, data_edges,groups=sample(c("A","B"),25,replace = TRUE)) %>%
          visNodes(borderWidth = 2) %>%
          visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
          visOptions (nodesIdSelection = TRUE,highlightNearest = TRUE)   %>%
          visPhysics(solver="forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant=-2000,centralGravity=0.02,springLength=40,springConstant=0.4,damping=1,avoidOverlap=1)) %>%
          visInteraction(multiselect = TRUE)#%>%
        
        
        
      })
      
      output$Cendroids<-renderText({
        if (is.null(mstValues$centroids)) {return()}
        mstValues$centroids
      })
      
      
      output$Links<-renderText({
        if (is.null(mstValues$keyvertices)) {return()}
        mstValues$keyvertices
      })
     
      
      
      observe({
        tempdata=list1()
        if (is.null(tempdata)) {return()}
        choice=list("Default")

        for (l in colnames(tempdata[[1]])){
          if (!is.numeric(tempdata[[1]][,l])){
            choice=c(choice,l)
            
          }}
        names(choice)=c("Default",unlist(choice)[-1])
        updateSelectInput(session,"colormst", label = "Select a Column for background", choices = choice)
        updateSelectInput(session,"bordermst", label = "Select a column for border", choices = choice)
      })
     
      
        
        observe({
          tempdata=(list1())
          x=input$colormst
          y=input$bordermst
          if (is.null(tempdata) | is.null(isolate(mstValues$flag2))){return()}
         
          if (x=="Default"){bgcolor=colors[mstValues$clusters]}
          else {bgcolor=colors[visualiseGenes(tempdata[[1]][,input$colormst])$label]}
          
          if (y=="Default"){bordcolor=colors[mstValues$clusters]}
          else  {bordcolor=colors[visualiseGenes(tempdata[[1]][,input$bordermst])$label]}
            
          visNetworkProxy("mstnetwork") %>%
              visUpdateNodes(data.frame(id=1:nrow(tempdata[[1]]),color.background=bgcolor,color.border=bordcolor))
          
        })
        
        centralValues=reactiveValues(centralities=data.frame(),rankings=data.frame,summary=data.frame())
        
        observeEvent(input$centralButton,{
            if (is.null(ig())){return()}
            igtemp=giant_component_extract(ig())[[1]]
            prop=proper_centralities(igtemp)
            central=calculate_centralities(igtemp,include=prop[c(-2,-27,-5,-6,-8,-21,-43,-38:-39,-31:-33)])
            central=central[lengths(central)!=0]
            x=as.data.frame(central,col.names = names(central),row.names=V(igtemp)$label)
           
            central_ranked=as.data.frame(apply(x,2,sort))
            colnames(central_ranked)=colnames(x)
          
            b=c("Average.Distance","Local.Bridging.Centrality","Wiener.Index.Centrality")
            
            #revert the order of indices 
            central_ranked[,-match(b,colnames(central_ranked))]=central_ranked[nrow(central_ranked):1,-match(b,colnames(central_ranked))]
          
            positions=c()
            for (i in 1:length(x))
            {positions_first=match(x[,i],central_ranked[,i])
            positions_last=match(x[,i],central_ranked[nrow(x):1,i])
            positions=cbind(positions,(positions_first-positions_last+nrow(x)+1)/2)
            }
            colnames(positions)=colnames(x)
            rownames(positions)=rownames(x)
            
            centralValues$centralities=x
            centralValues$rankings=positions
            centralValues$summary=as.data.frame(t(apply(positions,1,summary)),rownames=rownames(x))
            
            
          
        },ignoreInit = TRUE)
        
        output$Centralities<-DT::renderDataTable({
            if (is.null(centralValues$centralities)) {return()}
            centralValues$centralities
      
          
          })
       
         output$Rankings<-DT::renderDataTable({
            if (is.null(centralValues$rankings)) {return()}
            centralValues$rankings
          
          
         })
      
         
         output$Summary<-DT::renderDataTable({
           
           if (is.null(centralValues$summary)) {return()}
           centralValues$summary
           
          })
      
          
      
      })

    
        
       
       
       
   
          
        
        
        
      
      
         
       
         
         
         
         
         
       
          
      


  
  
  

