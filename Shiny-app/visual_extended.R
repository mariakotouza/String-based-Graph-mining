#confusionmatrix

confusion<-function(clustering1,clustering2,threshold=0.5){
  mat=table(clustering1,clustering2)
  mat1=prop.table(mat,1)
  mat2=t(prop.table(mat,2))
  x1=which(mat1>threshold,arr.ind = TRUE)
  x2=which(mat2>threshold,arr.ind = TRUE)
 
  
  x3=c()
  if (nrow(x1)>0 && nrow(x2)>0){
    for (i in 1:nrow(x1))
      {for (j in 1:nrow(x2))
        {
        temp=(x1[i,]==x2[j,c(2,1)])
             if (temp[1] && temp[2])
               x3=rbind(x3,x1[i,])
      }}}
   
  
  list("list1"=x1,"list2"=x2,"list3"=x3,"mat1"=mat1,"mat2"=mat2)
  }
