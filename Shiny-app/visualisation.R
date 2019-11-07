#visualise

visualiseGenes<-function(data)
  {udata=(unique(data))
  udata=udata[order(udata)]
  label=match(data,udata)
  list("name"=udata,"label"=label)}  