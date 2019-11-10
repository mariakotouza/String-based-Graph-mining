library(dplyr)

run_discretization <- function(datapath, dataset_name, numeric_file_name, ids_file_name, num_of_bins, num_of_bins_per_group,num_of_topics){
  data <- read.csv(paste0(datapath,"/", numeric_file_name),sep=" ",header = F)[,1:num_of_topics]

  if (file.exists(paste0(datapath,"/", ids_file_name))){
    data=cbind(sequence_id=read.csv(paste0(datapath,"/", ids_file_name),header = F),data)
  }else{
    data=cbind(sequence_id=1:nrow(data),data)
  }
  
  # normalization 
  data[,2:ncol(data)]=data[,2:ncol(data)]*10000/max(data[,2:ncol(data)]*10000)
  
  data_new=data
  data_discr=data
  
  #group values using ranges 0:1 with step=1/num_of_bins
  groups=seq(0,1,1/num_of_bins)
  #to do: make it more general ~ do not use letters that are only 24, use combinations
  let = LETTERS[1:num_of_bins] # The letters matrix
  let=let[length(let):1]
  sim=list()
  color=matrix(0,nrow=1, ncol=length(let))
  j=1
  overlaping <- F
  if (!overlaping){
    colors=colorRampPalette(c("blue", "green","yellow","red"))(num_of_bins/num_of_bins_per_group)
    for (i in 1:(num_of_bins/num_of_bins_per_group)){
      sim[[paste0(let[j],let[(j+num_of_bins_per_group-1)])]]=let[j:(j+num_of_bins_per_group-1)]
      color[j:(j+num_of_bins_per_group-1)]=colors[i]
      j=j+num_of_bins_per_group
    }
  }else{
    colors=colorRampPalette(c("blue", "green","yellow","red"))(num_of_bins-1)
    for (i in 1:(num_of_bins-1)){
      sim[[paste0(let[j],let[(j+num_of_bins_per_group-1)])]]=let[j:(j+num_of_bins_per_group-1)]
      color[j:(j+num_of_bins_per_group-1)]=colors[i]
      j=j+1
    }
  }
  
  num_of_topics <- ncol(data) - 1
  a <- list()
  max_group_length <- c()
  for (i in 1:num_of_topics){
    a[[i]] <- sim
    max_group_length <- c(max_group_length,length(a[[i]]))
  }
  
  sim <- a
  max_group_length <- max(max_group_length)
  
  #altsim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
  altsim <- sim
  # Naming the group of similarities
  for (i in 1:length(altsim)){
    names(altsim[[i]]) <- c("f","w","a","s","p","g","y","h","b","c","m")[1:(num_of_bins/num_of_bins_per_group)]
  }
  
  # Create custom colour scheme
  cs1 = make_col_scheme(chars=let,
                        cols=c(color))
  
  # Compute Letter Similarity 
  letter_sim=as.data.frame(matrix(0,nrow = length(let),ncol = length(let)))
  for (i in 1:length(let)){
    for (j in 1:length(let)){
      letter_sim[i,j]=1-abs(i-j)/length(let)
    }
  }
  colnames(letter_sim)=let
  row.names(letter_sim)=let
  
  for (i in 1:(length(groups)-1)){
    ids=which(data[,2:ncol(data)]>=groups[i] & data[,2:ncol(data)]<groups[i+1],arr.ind=TRUE)
    for (j in unique(ids[,2])){
      data_new[ids[which(ids[,2]==j),1],j+1]=let[i]
      data_discr[ids[which(ids[,2]==j),1],j+1]=groups[i]
    }
    if (i==length(groups)-1){
      ids=which(data[,2:ncol(data)]==groups[i+1],arr.ind=TRUE)
      for (j in unique(ids[,2])){
        data_new[ids[which(ids[,2]==j),1],j+1]= let[i] #let[i+1]
        data_discr[ids[which(ids[,2]==j),1],j+1]= groups[i]#let[i+1]
      }
    }
  }
  
  data_new$x <- apply( data_new[ , 2:ncol(data_new) ] , 1 , paste , collapse = "" )
  udata <- data_new[c(1,ncol(data_new))]
  colnames(udata)=c("Sequence.ID","AA.JUNCTION")
  
  tmp_path <- getwd()
  if (!file.exists(paste0(tmp_path,"/output_discr"))) dir.create(paste0(tmp_path,"/output_discr"))
  if (!file.exists(paste0("output_discr/", dataset_name))) dir.create(paste0("output_discr/", dataset_name))
  
  write.table(udata, file = paste0("output_discr/", dataset_name,"/udata.txt"), row.names = F, sep = "\t")
  write.table(x = read.csv(paste0(datapath,"/Topic_Similarity.txt"), sep="\t"), file = paste0("output_discr/", dataset_name,"/Topic_Similarity.txt"), sep = "\t", row.names = F)
  write.table(x = read.csv(paste0(datapath,"/true_classes.txt"), sep="\t"), file = paste0("output_discr/", dataset_name,"/true_classes.txt"), sep = "\t", row.names = F)
  
  udata_classes <- cbind(dataName=dataset_name, udata, Class=read.csv(paste0(datapath,"/true_classes.txt"), sep="\t", header = F), freq_cluster_id=1)
  write.table(x = udata_classes, file = paste0("output_discr/", dataset_name,"/udata_classes.txt"), sep = "\t", row.names = F)

  udata_classes_grouped <- udata_classes %>% group_by(dataName, AA.JUNCTION, V1) %>% dplyr::summarise(freq_cluster_id=n())
  udata_classes_grouped <- cbind(dataName=udata_classes_grouped$dataName, Sequence.ID=1:nrow(udata_classes_grouped),
                                 AA.JUNCTION=udata_classes_grouped$AA.JUNCTION,V1=udata_classes_grouped$V1,
                                 freq_cluster_id=udata_classes_grouped$freq_cluster_id)
  
  write.table(x = udata_classes_grouped, file = paste0("output_discr/", dataset_name,"/udata_classes_grouped.txt"), sep = "\t", row.names = F)
  
}