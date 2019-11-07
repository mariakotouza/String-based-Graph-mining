#clustering based on mst


mstClustering<-function(ig,threshold=1/15*dim(ig[])[1]){
  
  
  
  
  
  
   mstig=mst(ig,algorithm='prim')
  # dis=distances(mstig)
   degree=centr_degree(mstig)
   centroids=which(degree$res>threshold)
  # 
  # 
  # 
  # mins=apply(dis[,centroids],1,which.min)
  # clusters=match(centroids[mins],centroids)
  # 
  # 
  # df=as_data_frame(mstig,what = "edges")
  # linkedges=which(clusters[df[,"from"]]!=clusters[df[,"to"]])
  # keyvertices=df[linkedges,c("from","to")]
  # 
  # 
  # list("clusters"=clusters,"centroids"=centroids,"keyvertices"=keyvertices)
   return (centroids)
}
