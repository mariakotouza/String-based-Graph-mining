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
    imp_node_ids <- which(c_s > median(c_s))
    for (i in imp_node_ids){
      # for rows
      j <- which(w[i,]!=0)
      w[i,j] <- w[i,j] + c_s[i]
      
      # for cols
      j <- which(w[,i]!=0)
      w[j,i] <- w[j,i] + c_s[i]
      print(c_s[i])
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
