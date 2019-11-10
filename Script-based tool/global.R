
#filterg=data.frame("Columns"=character(),"Keys"=character(),"I/E"=character(),stringsAsFactors = FALSE)
#colors=colors()[c(85:131)*5]
library(dplyr)
library(shiny) 
library(stringdist)
library(igraph)
library(visNetwork)
library(CINNA)
library(DT)
library(ggraph)
library(graphlayouts)
library(cluster)
library(shinyalert)
library(shinybusy)
library(shinythemes)
library(optrees)
library(shallot)
library(aricode)
library(gplots)
library(shinyFiles)
library(shinyBS)
library(rfUtilities)


source('stringdistances.R', echo=TRUE)
source('visualisation.R', echo=TRUE)
source('filtered.R', echo=TRUE)
source('multiple_clustering.R', echo=TRUE)
source('visual_extended.R', echo=TRUE)
source('clusterMetrics.R', echo=TRUE)
source('MST.R')

tmp_path<-getwd() #change it to "/tmp" for server
#logfile
if(!file.exists(paste0(tmp_path,"/log_files"))){ 
  dir.create(paste0(tmp_path,"/log_files"))
}
logFile = paste0(tmp_path,"/log_files/log_file ",trunc(as.numeric(Sys.time())),".txt")
cat(paste0("Function","\t","Parameters","\t","Num of input rows","\t","Num of input columns","\t","Start time","\t","End time","\t","Memory used"), file=logFile, append=FALSE, sep = "\n")

save_tables_individually=T

colors=unique(gsub('[0-9]+','', colors()))
colors=colors[length(colors):1]
# colors=c("green","yellow","red","blue","orange","purple","grey","pink","black","brown","lightblue")
# colors=rep(colors,19)
id=NULL

domain <- "AA"  #items
if (domain == "AA"){
  sequences_for_distance_calc <- list("IMGT.gapped.nt.sequences.V.D.J.REGION"="IMGT.gapped.nt.sequences.V.D.J.REGION",
                                      "IMGT.gapped.AA.sequences.V.D.J.REGION"= "IMGT.gapped.AA.sequences.V.D.J.REGION",
                                      "IMGT.gapped.nt.sequences.V.J.REGION"="IMGT.gapped.nt.sequences.V.J.REGION", 
                                      "IMGT.gapped.AA.sequences.V.J.REGION"="IMGT.gapped.AA.sequences.V.J.REGION",
                                      "Summary.AA.JUNCTION"="Summary.AA.JUNCTION")
  default_seq_for_distance_calc <- "Summary.AA.JUNCTION"
}else{
  sequences_for_distance_calc <- list("AA.JUNCTION"="AA.JUNCTION")
  default_seq_for_distance_calc <- "AA.JUNCTION"
}

cols_data <- c()