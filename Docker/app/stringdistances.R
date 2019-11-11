
stringdistances<-function(seq,algo){
  sim=stringdistmatrix(seq,method=algo,useBytes = FALSE,q=2)
  if (algo %in% c("dl","lv","osa","lcs","hamming","qgram")){
    if (algo %in% c("lcs","qgram"))
      funn=sum
    else 
      funn=max
    lengths=nchar(as.character(seq))
    normal=as.vector(combn(lengths,2,FUN=funn))
    sim=sim/normal
  }
  return (sim)
}