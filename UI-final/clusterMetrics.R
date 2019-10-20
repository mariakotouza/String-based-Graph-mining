#clustering metrics 



conductance<-function(graph,membership)
{
  edges=as_data_frame(graph)
  edges$from=membership[edges$from]
  edges$to=membership[edges$to]
  con=c()
  total=sum(edges$weight)
  
  conduct<-function(i,edges,total)
  {     
    inter=sum(edges$weight[xor(edges$to==i,edges$from==i)])
    clust=sum(edges$weight[edges$to==i | edges$from==i])
    con=inter/min(total-clust,clust)
    return(con)
  }
  
  con<-lapply(sort(unique(membership)),conduct,edges=edges,total=total)
  
  
  list (conductance=1-mean(unlist(con)),conductances=con)
}




coverage<-function(graph,membership)
{
  edges=as_data_frame(graph)
  edges$from=membership[edges$from]
  edges$to=membership[edges$to]
  return(sum(edges$weight[edges$from==edges$to])/sum(edges$weight))
}
