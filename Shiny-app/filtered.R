#filtering processes

#install.packages("dplyr")
#library("dplyr")


uniteDatasets<-function(...)
{
  x=list(...)
  dataset=c()
  for (i in 1:length(x))
    dataset=rbind(dataset,x[i])
  dataset=as.data.frame(unique(dataset))
  return (dataset)
}

subtractDatasets<-function(minuend,subtracter)
{
  minuend=setdiff(minuend,subtracter)
  return (minuend)
}


# includeInGraph<-function(dataset,Summary.D.GENE.and.allele=NULL,Summary.J.GENE.and.allele=NULL,Summary.V.GENE.and.allele=NULL)
#   {
#   data=dataset
#   if (!is.null(Summary.D.GENE.and.allele)){
#       data=c()
#       for (i in 1:length(Summary.D.GENE.and.allele))
#         data=rbind(data,dataset[dataset$Summary.D.GENE.and.allele==Summary.D.GENE.and.allele[i],])
#   }    
#   
#   dataset=data
#   
#   if (!is.null(Summary.J.GENE.and.allele)){
#     data=c()
#     for (i in 1:length(Summary.J.GENE.and.allele))
#       data=rbind(data,dataset[dataset$Summary.J.GENE.and.allele==Summary.J.GENE.and.allele[i],])
#   }   
#   
#   if (!is.null(Summary.V.GENE.and.allele)){
#     data=c()
#     for (i in 1:length(Summary.V.GENE.and.allele))
#       data=rbind(data,dataset[dataset$Summary.V.GENE.and.allele==Summary.V.GENE.and.allele[i],])
#   } 

# return(data)
#  }



includeInGraph<-function(data,x)
{
  
  data=cbind(data,id=c(1:nrow(data)))
  
  for (l in x)
  {
    tempdata=c()
    if (!is.numeric(data[,l[1]]))
      tempdata=data[data[,l[1]] %in% l[2:length(l)],]
    
    else 
      tempdata=data[(data[,l[1]]>l[2] & data[,l[1]]<l[3] ),]
    data=tempdata  
  }
  return (sort(data$id))   
}





excludeFromGraph<-function(data,x)
{
  
  data=cbind(data,id=c(1:nrow(data)))
  
  for (l in x)
    if (!is.numeric(data[,l[1]]))
      data=data[-which(data[,l[1]] %in% l[2:length(l)]),]
    else
      data=data[-which(data[,l[1]]>l[2] & data[,l[1]]<l[3]),]
    return (data$id)
}



