# Evaluation metrics -> depends on the domain
evaluation_metrics <- c("ff", "Num_of_clusters", "Average-Identity", "Average-Similarity", "Average-Entropy_Sim",
                        "Average-BS", "Average-TopicSim","SD-Similarity", "SD-Identity","Average-Entropy_Id", 
                        "SD-Entropy_Id", "SD-Entropy_Sim", "SD-BS", "SD-TopicSim")

flagtic <- T

random_data_generator <- function(n,l=20){
  a <- do.call(paste0, replicate(l, sample(let, n, TRUE), FALSE))
  a
}

#data_let: data frame with columns ID, AA.JUNCTION
let_to_arithm <- function(data_let, num_of_items){
  gen_data=data.frame(ID=1:num_of_items,AA.JUNCTION=data_let,stringsAsFactors = F)
  gen_data$ID=row.names(gen_data)
  
  a=as.data.frame(str_split_fixed(gen_data$AA.JUNCTION, "", num_of_topics),stringsAsFactors = F)
  for (i in 1:ncol(a)){
    a[,i]<-factor(a[,i],levels=let)
  }
  a=ddply(a,colnames(a),function(x){as.numeric(x)/num_of_bins-1/num_of_bins})
  data_arithm=cbind(ID=1:num_of_items,a)
  data_arithm
}

####################### Binary Tree Construction Functions #########################
#list1 = listb
#list1 = result2
# Function Matrices compute the apsolute matrix, the matrix with percentages, the absolute similarity matrix and the percentage similarity matrix
# Function Matrices contains the old function Finish. Finish compute cluster's identity and similarity and stop cluster's division if identity or similarity is greater than endper
Matrices <- function(list1,let,sim,d,algo,algocol,backcol,backcolj6,logFile,domain,max_group_length){
  tic()
  udata = list1$udata
  br = list1$br
  permat = list1$permat
  persim = list1$persim
  cl = list1$cl
  listxx = list1$list
  listyy = list1$listn
  sumper = list1$sumper # The total percentage of permat 
  sumper2 = list1$sumper2 # The total percentage of persim
  endper = list1$endper
  dfsum = list1$dfsum
  nn = list1$nn
  last = list1$last
  progend  = list1$progend
  leaf = list1$leaf
  leaf2 = list1$leaf2
  
  if (domain!="AA"){
    mymat = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
    permat = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
    simmat = matrix(0,nrow = max_group_length,ncol =str_length(udata$AA.JUNCTION[1]))
    persim = matrix(0,nrow = max_group_length+1,ncol =str_length(udata$AA.JUNCTION[1]))
    rownames(mymat) = let
    rownames(permat) = c(let,"Entropy")
    rownames(simmat) = c(paste0("group_",1:max_group_length))
    rownames(persim) = c(paste0("group_",1:max_group_length),"Entropy")
  }else{
    # Find the sequences with gene J6
    ind1 = str_which(udata[udata$clusters == br,]$J.GENE.and.allele,"J6")
    qw = 1:length(udata[udata$clusters == br,]$AA.JUNCTION)
    if(length(ind1) == 0){
      ind2 = qw
    }else{
      ind2 = qw[-ind1]
    }
    
    # Initialize matrices for J6 sequences
    mymat1 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
    simmat1 = matrix(0,nrow = max_group_length,ncol =str_length(udata$AA.JUNCTION[1]))
    rownames(mymat1) = let
    rownames(simmat1) = c(paste0("group_",1:max_group_length))
    
    # Initialize matrices for NO J6 sequences
    mymat2 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
    simmat2 = matrix(0,nrow = max_group_length,ncol =str_length(udata$AA.JUNCTION[1]))
    rownames(mymat2) = let
    rownames(simmat2) = c(paste0("group_",1:max_group_length))
    
    # Initialize matrices for all sequences
    permat = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
    rownames(permat) = c(let,"Entropy")
    persim = matrix(0,nrow = max_group_length+1,ncol =str_length(udata$AA.JUNCTION[1]))
    rownames(persim) = c(paste0("group_",1:max_group_length),"Entropy")
    
    # Compute matrices for J6 sequences
    if(length(ind1) != 0 ){
      trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION[ind1],"")
      align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
      for(i in (algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)){
        temptab = plyr::count(align[i],vars = colnames(align)[i])
        names(temptab)[1] = "X1"
        mymat1[which(is.na(match(names(mymat1[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
        temptab1 = temptab
        gg = as.vector(unlist(temptab1[1]))
        temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),gg[x]))])))
        #temptab1[1] = names(sim[[i]])[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim,x)))]
        temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
        simmat1[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),i] = as.vector(unlist(temptab1[2]))
      }
    }
    
    # Compute matrices for NO J6 sequences
    if(length(ind2) != 0 ){
      trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION[ind2],"")
      align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
      for(i in (algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)){
        temptab = plyr::count(align[i],vars = colnames(align)[i])
        names(temptab)[1] = "X1"
        mymat2[which(is.na(match(names(mymat2[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
        permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
        temptab1 = temptab
        gg = as.vector(unlist(temptab1[1]))
        temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim[[i]])[which(str_detect(unlist(sim[[i]]),gg[x]))]))) #delete digits from the names of similarity groups
        temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
        simmat2[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),i] = as.vector(unlist(temptab1[2]))
        persim[max_group_length+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
      }
    }
    
    # Combine matrices for No J6 and J6 sequences
    mymat = mymat1 + mymat2
    simmat = simmat1 + simmat2
    meg = length(ind1) + length(ind2)
    if(length(ind2) == 0 || backcol == backcolj6){
      permat[1:length(let),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (mymat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
      persim[1:max_group_length,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (simmat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
    }else{
      permat[1:length(let),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (mymat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
      persim[1:max_group_length,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (simmat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
      if(backcolj6 != 0){
        permat[1:length(let),(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] = (mymat[,(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] / length(ind2)) * 100
        persim[1:max_group_length,(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] = (simmat[,(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] / length(ind2)) * 100
      }
    }
  }
  
  # Find the Entropy for all the sequences
  trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION,"")
  align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
  
  for(i in 1:str_length(udata[udata$clusters == br,]$AA.JUNCTION[1])){
    temptab = plyr::count(align[i],vars = colnames(align)[i])
    names(temptab)[1] = "X1"
    match_elements=match(names(mymat[,i]),as.vector(unlist(temptab[1])))
    ids=as.vector(unlist(temptab[1]))[na.omit(match_elements)]
    mymat[ids,i] = as.vector(unlist(temptab[2]))[na.omit(match_elements)]
    #mymat[which(is.na(match(names(mymat[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
    permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
    temptab1 = temptab
    temptab1[1] = names(sim[[i]])[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim[[i]],x)))]
    temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
    
    simmat[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),i] = as.vector(unlist(temptab1[2]))
    persim[max_group_length+1,i] = entropy(simmat[,i],base = exp(1))
  }
  permat[1:length(let),] = (mymat / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100
  persim[1:max_group_length,] = (simmat / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100
  
  if (domain=="AA"){
    if(algocol != 0 && backcol!=0){
      # Initialize matrices
      mymat3 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
      simmat3 = matrix(0,nrow = max_group_length,ncol =str_length(udata$AA.JUNCTION[1]))
      rownames(mymat3) = let
      rownames(simmat3) = c(paste0("group_",1:max_group_length))
      permat3 = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
      rownames(permat3) = c(let,"Entropy")
      persim3 = matrix(0,nrow = max_group_length+1,ncol =str_length(udata$AA.JUNCTION[1]))
      rownames(persim3) = c(paste0("group_",1:max_group_length),"Entropy")
      
      trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION,"")
      align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
      exc = 1:algocol
      exc2 = (str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol):str_length(udata[udata$clusters == br,]$AA.JUNCTION[1])
      exc = c(exc,exc2)
      for(i in 1:length(exc)){
        temptab = plyr::count(align[exc[i]],vars = colnames(align)[exc[i]])
        names(temptab)[1] = "X1"
        mymat3[which(is.na(match(names(mymat3[,exc[i]]),as.vector(unlist(temptab[1])))) == FALSE),exc[i]] = as.vector(unlist(temptab[2]))
        permat3[length(let) + 1,exc[i]] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
        temptab1 = temptab
        gg = as.vector(unlist(temptab1[1]))
        temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),gg[x]))])))
        #temptab1[1] = names(sim[[i]])[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim,x)))]
        temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
        simmat3[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),exc[i]] = as.vector(unlist(temptab1[2]))
        persim3[max_group_length+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
      }
      permat3[1:length(let),] = (mymat3 / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100
      persim3[1:max_group_length,] = (simmat3 / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100
      permat[1:length(let),exc] = permat3[1:length(let),exc] 
      persim[1:max_group_length,exc] = persim3[1:max_group_length,exc] 
    }
  }
  
  listxx$temp = permat
  names(listxx)[length(listxx)] = sprintf('permat_br.%d', br) # Save the permat with this format
  listyy$temp = persim
  names(listyy)[length(listyy)] = sprintf('persim_br.%d', br)
  
  ################################ Finish ###########################
  t1 =which(permat[-length(permat),] == 100,arr.ind = TRUE)
  sumper = (length(as.numeric(t1[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  #find the number of non-zero groups. if this number is = 2 than sim=100
  t2 =which(persim[-length(persim),] >= 100,arr.ind = TRUE)
  
  sumper2 = (length(as.numeric(t2[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  vv = length(udata[udata$clusters == br,]$AA.JUNCTION)
  dfsum[nrow(dfsum) + 1,] = c(sumper,sumper2,br,vv)
  if(algo == "Identity"){
    if (sumper > endper && leaf == FALSE){
      nn = TRUE  # When nn = TRUE the percentage of sumper < endper%
      if (sumper2 > endper && leaf == FALSE){
        leaf2 = TRUE			
      }
    }
  }else{
    if (sumper2 > endper){
      leaf2 = TRUE	
    }
  }
  
  result1 = list("listax" = list1$listax,"ggdf" = list1$ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep"= list1$ep, "clep" = list1$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "selected_cell_id" = list1$selected_cell_id, "cel" = list1$cel,"endper" = list1$endper, "last" = last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  
  en = toc(quiet = TRUE)
  cat(paste0("Matrices","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  if ( leaf == TRUE){
    nn = FALSE
    result1 = list("listax" = list1$listax,"ggdf" = list1$ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep"= list1$ep, "clep" = list1$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "selected_cell_id" = list1$selected_cell_id, "cel" = list1$cel,"endper" = list1$endper, "last" = last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
    return(result1)
  }else{
    return(result1)
  }
}


#list2 = result1
# Function Choice chooses which matrix cell will be used for the division of the data 
# Function Choice contains the old Divide. Divide is responsible for dividing the sequences based on the cell given.
# Function Choice also contains the old Control. Control is responsible for updating level, cluster and other counters as well as for stopping the algorithm 
Choice <- function(list2,let,sim,d,algo,algocol,backcol,backcolj6,logFile,domain,max_group_length){
  # Count the execution time of Choices
  tic()
  # Keep data up to date
  udata = list2$udata
  br = list2$br
  cl = list2$cl
  nn = list2$nn
  cel = list2$cel
  selected_cell_id = list2$selected_cell_id
  met = list2$met
  ep = list2$ep
  clep = list2$clep
  ggdf = list2$ggdf
  progend  = list2$progend
  leaf = list2$leaf
  listax = list2$listax
  leaf2 = list2$leaf2
  
  # Find desired cell for division
  if(nn == TRUE || (algo == "Similarity") ){
    permat = list2$persim # If sumper < endper% we want to check only the persim matrix
  }else{
    permat =list2$permat  # Else permat and if it is necessary the persim matrix
    persim = list2$persim
  }
  
  if (domain!="AA"){
    permat_temp <- permat
  }else{
    permat_temp <- permat[,(algocol + 1):(ncol(permat)-backcol)]
  }
  
  cel = which(permat_temp == max(permat_temp), arr.ind = TRUE)
  selected_cell_id = 1
  poss = max(permat_temp) 
  # We exclude the 100 % from the max values
  if (max(permat_temp) == 100){ 
    cel = which(permat_temp == max(permat_temp[permat_temp!=max(permat_temp)]), arr.ind = TRUE) # The desired cell
    poss = max(permat_temp[permat_temp!=max(permat_temp)])
  }
  
  # Check if clusters division is stopped or not
  if(leaf2 == FALSE && poss != 0){  
    # if cel contain more than one cells, find the best cell matching some criteria
    selected_cell_id = 1
    if(algo == "Identity"){
      if ((length(cel)/2) > 1){
        dddff = min(permat_temp[nrow(permat_temp),cel[,2]])
        dddff2 = which(permat_temp[nrow(permat_temp),cel[,2]] == dddff)
        selected_cell_id = dddff2[1]
        if (length(dddff2)>1 && nn == FALSE){ # If the vector has 2 or more numbers means that we have columns with the same entropy and nn = FALSE in order not to double check the persim
          dddff3 = persim[,(algocol + 1):(ncol(permat)-backcol)][sapply(1:length(dddff2), function (x){ str_which(names(sim[[selected_cell_id]]),gsub('[[:digit:]]+', '', names(unlist(sim[[selected_cell_id]])[which(str_detect(unlist(sim[[selected_cell_id]]),let[cel[x,1]]))])))}),cel[,2]]
          dddff3 = max(diag(dddff3))
          dddff4 = which(diag( persim[,(algocol + 1):(ncol(permat)-backcol)][sapply(1:length(dddff2), function (x){str_which(names(sim[[selected_cell_id]]),gsub('[[:digit:]]+', '', names(unlist(sim[[selected_cell_id]])[which(str_detect(unlist(sim[[selected_cell_id]]),let[cel[x,1]]))])))}),cel[,2]]) == dddff3)
          selected_cell_id = dddff4[1]
          if(length(dddff4)> 1){
            dddff5 = min(persim[,(algocol + 1):(ncol(permat)-backcol)][nrow(persim[,(algocol + 1):(ncol(permat)-backcol)]),cel[dddff4,2]])
            dddff6 = which(persim[,(algocol + 1):(ncol(permat)-backcol)][nrow(persim[,(algocol + 1):(ncol(permat)-backcol)]),cel[dddff4,2]] == dddff5)
            selected_cell_id = dddff6[1] 
          }
        }
      }
    }else{
      if ((length(cel)/2) > 1){
        dddff = min(permat_temp[nrow(permat_temp),cel[,2]])
        dddff2 = which(permat_temp[nrow(permat_temp),cel[,2]] == dddff)
        selected_cell_id = dddff2[1]
      }
    }
    nn = FALSE # Return nn in it's original value 
    
    ##################################################### Old Divide ############################
    # If we need a new level, then we create a new column in udata dataframe with its name (level.ep)
    if (met == 0) {
      udata$temp = NA
      names(udata)[length(udata)] = sprintf('level.%d', ep)
    }
    
    if (algo == "Identity"){
      # Find sequences that contain the cell's letter in cell's position 
      x1 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), let[cel[selected_cell_id,1]]), "TRUE")
      mk1 = length(x1)
      cltp1 = cl + 1
      y1 = udata[udata$clusters == br,]$AA.JUNCTION
      z1 = y1[x1]
      
      # The other sequences of the cluster
      x2 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), let[cel[selected_cell_id,1]]), "FALSE")
      mk2 = length(x2)
      cltp2 = cl + 2
      y2 = udata[udata$clusters == br,]$AA.JUNCTION
      z2 = y2[x2]
      lengdif = FALSE
      
      # Find wich new sub-cluster has more sequences and give it first cluster name
      if(mk1 < mk2){
        lengdif = TRUE
        templeng = z1
        z1 = z2
        z2 = templeng
        temp2 = x1
        x1 = x2
        x2 = temp2
        ggdf[nrow(ggdf) + 1,] = c(cltp1,mk2)
        ggdf[nrow(ggdf) + 1,] = c(cltp2,mk1)
      }else{
        ggdf[nrow(ggdf) + 1,] = c(cltp1,mk1)
        ggdf[nrow(ggdf) + 1,] = c(cltp2,mk2)
      }
      
      listax$temp = as.vector(unlist(listax[sprintf('cl.%d', br)]))[x1]
      names(listax)[length(listax)] = sprintf('cl.%d', cl+1) 
      listax$temp = as.vector(unlist(listax[sprintf('cl.%d', br)]))[x2]
      names(listax)[length(listax)] = sprintf('cl.%d', cl+2) 
      udata[as.vector(unlist(listax[sprintf('cl.%d', br)]))[x1],sprintf('level.%d', ep)] = cl+1
      udata[as.vector(unlist(listax[sprintf('cl.%d', br)]))[x2],sprintf('level.%d', ep)] = cl+2
      
      if(lengdif == TRUE){
        udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), let[cel[selected_cell_id,1]]), cl+2 ,cl+1)
      }else{
        udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), let[cel[selected_cell_id,1]]), cl+1 ,cl+2)
      }
    }else{
      # Find sequences contains the cell's similarity group in cell's position 
      strings.to.find = unlist(sim[[selected_cell_id]][cel[selected_cell_id,1]])
      x1 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), str_c(strings.to.find, collapse="|")), "TRUE")
      mk1 = length(x1)
      cltp1 = cl + 1
      y1 = udata[udata$clusters == br,]$AA.JUNCTION
      z1 = y1[x1]
      
      # The other sequences of the cluster
      x2 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), str_c(strings.to.find, collapse="|")), "FALSE")
      mk2 = length(x2)
      cltp2 = cl + 2
      y2 = udata[udata$clusters == br,]$AA.JUNCTION
      z2 = y2[x2]
      lengdif = FALSE
      
      # Find wich new sub-cluster has more sequences and give it first cluster name
      if(mk1 < mk2){
        lengdif = TRUE
        templeng = z1
        z1 = z2
        z2 = templeng
        temp2 = x1
        x1 = x2
        x2 = temp2
        ggdf[nrow(ggdf) + 1,] = c(cltp1,mk2)
        ggdf[nrow(ggdf) + 1,] = c(cltp2,mk1)
      }else{
        ggdf[nrow(ggdf) + 1,] = c(cltp1,mk1)
        ggdf[nrow(ggdf) + 1,] = c(cltp2,mk2)
      }
      
      listax$temp = as.vector(unlist(listax[sprintf('cl.%d', br)]))[x1]
      names(listax)[length(listax)] = sprintf('cl.%d', cl+1) # Save the permat with this format
      listax$temp = as.vector(unlist(listax[sprintf('cl.%d', br)]))[x2]
      names(listax)[length(listax)] = sprintf('cl.%d', cl+2) # Save the permat with this format
      udata[as.vector(unlist(listax[sprintf('cl.%d', br)]))[x1],sprintf('level.%d', ep)] = cl+1
      udata[as.vector(unlist(listax[sprintf('cl.%d', br)]))[x2],sprintf('level.%d', ep)] = cl+2
      
      if(lengdif == TRUE){
        udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), str_c(strings.to.find, collapse="|")), cl+2 ,cl+1) 
      }else{
        udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), str_c(strings.to.find, collapse="|")), cl+1 ,cl+2) 
      }
    }
    clep[cl+1] = ep # Level of the cluster cl+1  
    clep[cl+2] = ep # Level of the cluster cl+2
    cl = cl + 2 # Increase the cluster by 2
  }else{
    leaf2 = FALSE
    nn = FALSE
  }
  
  ############################################ Old Control #############################################
  br = br + 1 # Increase the branch by 1
  met = met + 2 # Increase the counter by 2
  list2 = list("listax" = listax,"ggdf" = ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "selected_cell_id" = selected_cell_id, "cel" = cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  if( ((clep[br-1] < clep[br]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE && (length(which(udata$clusters == br)) > 2) ){ # If the next branch is in the next level
    met = geomSeq(1,2,1,1000)[ep+1]
  }
  while (length(which(udata$clusters == br)) <= 1){ # While the number of sequences in the branch is less than 1, go to the next branch and change counter 
    # Condition wich ends the algorithm 
    if(is.na(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]))){
      clmax = cl
      brtemp = br
      progend = TRUE
      list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "selected_cell_id" = list2$selected_cell_id, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
      if(brtemp < clmax){
        # Find the statistics of the remaining clusters before you end
        for (i in (br+1):clmax) {
          brtemp = i
          list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "selected_cell_id" = list2$selected_cell_id, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
          list2 = Matrices(list2,let,sim,d,algo,algocol,backcol,backcolj6,logFile,domain,max_group_length)
        }
      }
      en = toc(quiet = TRUE)
      cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
      return(list2)
    }
    
    # Run Matrices as leaf
    leaf = TRUE
    list2 = list("listax" = listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "selected_cell_id" = selected_cell_id, "cel" = cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
    list2 = Matrices(list2,let,sim,d,algo,algocol,backcol,backcolj6,logFile,domain,max_group_length)
    
    # Condition wich ends the algorithm 
    if(is.na(((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))))){
      clmax = cl
      brtemp = br
      progend = TRUE
      list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "selected_cell_id" = list2$selected_cell_id, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
      if(brtemp < clmax){
        # Find the statistics of the remaining clusters before you end
        for (i in (br+1):clmax) {
          brtemp = i
          list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "selected_cell_id" = list2$selected_cell_id, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
          list2 = Matrices(list2,let,sim,d,algo,algocol,backcol,backcolj6,logFile,domain,max_group_length)
        }
      }
      en = toc(quiet = TRUE)
      cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
      return(list2)
    }
    
    # If the next branch is in the next level
    if( ((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE ){
      met = 0
      ep = ep + 1
      br = br + 1
    }else{
      br = br + 1
      met = met +2 
    }
    list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "selected_cell_id" = list2$selected_cell_id, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  }
  
  # When the counter reaches the end value (geometric sequence) we increase the level counter
  if( met == geomSeq(1,2,1,1000)[ep+1]){ 
    met = 0
    ep = ep + 1
  }
  
  # The return list
  leaf = FALSE
  result2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "selected_cell_id" = selected_cell_id, "cel" = cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  en = toc(quiet = TRUE)
  cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
  return(result2)
  
}


# A function that generates a geometric sequence
geomSeq <- function(start,ratio,begin,end){
  begin=begin-1
  end=end-1
  start*ratio**(begin:end)
}

########################  Other Functions ######################## 
# Creating a Static Tree visualization
Den <- function(lev,df,flagtic,logFile){
  if(flagtic == TRUE) tic(sprintf("Den -- Level: %d ", lev))
  df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
  df_args <- c(df[str_which(names(df), "level.")], sep="/")
  if(lev == max(na.omit(clep))){
    df$pathString<- do.call(paste, df_args)
    kk = df$pathString
    for(i in 1:length(kk)){
      temp = str_locate(kk[i],"/NA")
      if(is.na(temp[1]) == FALSE){
        temp2 = str_sub(kk[i], 1, temp[1]-1);
        kk[i] = temp2 
      }
    } 
  }else{
    df$pathString<- do.call(paste, df_args)
    kk = df$pathString
    gg =as.data.frame(str_locate_all(kk,"/"))
    tem = 1
    for (i in 1:length(kk)){
      kk[i] = str_sub(kk[i],1,gg[,tem][lev+1]-1)
      tem = tem + 2
    }
    for(i in 1:length(kk)){
      temp = str_locate(kk[i],"/NA")
      if(is.na(temp[1]) == FALSE){
        temp2 = str_sub(kk[i], 1, temp[1]-1);
        kk[i] = temp2 
      }
    } 
  }
  df$pathString = kk
  x <- ToDataFrameTree(df, "pathstring")
  #par = str_count(df$pathString,"/")
  #str_locate(df$pathString,"/")
  #par = rev(gregexpr("\\/", df$pathString))
  xN <- as.Node(x)
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  plot(xN)
}

# A function to plot logos of level and leaf (until this level) clusters
LogoLev <- function(lev,Clus,df,flagtic,logFile){
  if(flagtic == TRUE) tic(sprintf("LogoLev -- Level: %d ", lev))
  if( is.element(lev,clep) == FALSE){
    print("Den yparxei")
  }else{
    t1 = which(clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    listff = list()
    if(length(t1) %% 3 == 0){
      nc = length(t1)%/%3
    }else{
      nc = length(t1)%/%3 + 1
    }
    for(i in 1:length(t1)){
      if(i<= length(xm2)){
        listff$temp = na.omit(df[df[names(df) == sprintf('level.%d', clep[t1[i]])] == t1[i],]$AA.JUNCTION)
        names(listff)[length(listff)] = sprintf('Cluster.%d - num:%d - leaf', t1[i],Clus$seqnum[t1[i]+1]) 
      }else{
        x1 = as.data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lev))] == t1[i], ]$AA.JUNCTION))
        if (nrow(x1) >0){
          names(x1)[1]= "AA.JUNCTION"
          x1 = as.character(x1$AA.JUNCTION)
          listff$temp = x1
          names(listff)[length(listff)] = sprintf('Cluster.%d - num:%d', t1[i],Clus$seqnum[t1[i]+1]) 
        } 
      }
    }
    
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggseqlogo(listff, ncol=nc, method = "prob",col_scheme=cs1) 
  }
} 

LogoLevNew <- function(lev){
  if(flagtic == TRUE) tic(sprintf("LogoLev -- Level: %d ", lev))
  if( is.element(lev,clep) == FALSE){
    print("Den yparxei")
  }else{
    t1 = which(clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    
    if(seqen == TRUE){
      if(length(which(Clus[t1+1,]$seqnum < seqthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$seqnum < seqthr)]
      }
    }else if(ideen == TRUE){
      if(length(which(Clus[t1+1,]$Identity < idethr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Identity < idethr)]
      }
    }else if (simen == TRUE){
      if(length(which(Clus[t1+1,]$Similarity < simthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Similarity <  simthr)]
      }
    }else if (allen == TRUE){
      if(length(which(Clus[t1+1,]$seqnum < seqthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$seqnum < seqthr)]
      }
      if(length(which(Clus[t1+1,]$Identity < idethr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Identity < idethr)]
      }
      if(length(which(Clus[t1+1,]$Similarity <  simthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Similarity <  simthr)]
      }
    }
    
    #swap(a$x,a$freq)
    #a = plyr::count(Clus[t1,]$Similarity)
    #barplot(t(as.matrix(a)), width=2, main = sprintf('Cluster.%d',cl)) 
    #legend("topright",inset=c(-0.03,0), fill=heat.colors(nrow(a)), legend=let,cex = 0.6)
    #if(length(t1) %% 3 == 0){
    #  nc = length(t1)%/%3
    #}else{
    #  nc = length(t1)%/%3 + 1
    #}
    listff = list()
    nc = 3
    for(i in 1:length(t1)){
      if(t1[i]<= max(xm2)){
        listff$temp = na.omit(df[df[names(df) == sprintf('level.%d', clep[t1[i]])] == t1[i],]$AA.JUNCTION)
        names(listff)[length(listff)] = sprintf('Cluster.%d - num:%d - leaf', t1[i],Clus$seqnum[t1[i]+1]) 
      }else{
        x1 = as.data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lev))] == t1[i], ]$AA.JUNCTION))
        if (nrow(x1) >0){
          names(x1)[1]= "AA.JUNCTION"
          x1 = as.character(x1$AA.JUNCTION)
          listff$temp = x1
          names(listff)[length(listff)] = sprintf('Cluster.%d - num:%d', t1[i],Clus$seqnum[t1[i]+1]) 
        } 
      }
    }
    
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggseqlogo(listff, ncol=nc, method = "prob",col_scheme=cs1) 
  }
}

# A function to plot the logo of a cluster
LogoCl <- function(cl,df,flagtic,logFile){
  if(flagtic == TRUE) tic(sprintf("LogoCl -- Cluster: %d ", cl))
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  if(is.na(clep[cl])){
    print("Den yparxei")
  }else{
    ggseqlogo(na.omit(df[df[names(df) == sprintf('level.%d', clep[cl])] == cl,]$AA.JUNCTION), method = "prob", col_scheme=cs1)
  }
}

# Copmute the sequnces, Identity, Similarity of level and leaf (until this level) clusters and make barplots
SatLev <- function(lev,Clus,df,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if( is.element(lev,clep) == FALSE){
    NULL
  }else{
    xN <- as.Node(df)
    t1 = which(clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    par(mfrow=c(1,3))
    a = table(Clus[t1+1,]$seqnum)
    barplot(a, width=2, main = sprintf('Sequences for level.%d' ,lev))
    b = table(round(Clus[t1+1,]$Identity,2))
    barplot(b, width=2, main = sprintf('Identity of sequences for level.%d' ,lev))
    c = table(round(Clus[t1+1,]$Similarity,2))
    barplot(c, width=2, main = sprintf('Similarity of sequences for level.%d' ,lev))
  }
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0(sprintf("LogoLev -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
}
#df = x
# Creating a plot for cluster 5 for example

# A function to plot barplots of level and leaf (until this level) clusters
BarLev <- function(lev,clep,xN,perlist,Clus,let,flagtic,logFile,jpeg_label,output_folder,slash_for_topics,data_topic,algo,num_of_topics){
  if(flagtic == TRUE) tic(sprintf("BarLev -- Level: %d ", lev))
  if( is.element(lev,clep) == FALSE){
    print("Den yparxei")
  }else{
    jpeg(paste0(output_folder,slash_for_topics,data_topic,"/","gg_bar_plots_level_",lev,"_algo_",algo,jpeg_label,".jpg"), width=3900, height=6000, res=100)
    t2 = which(clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(clep == lev))
    xm2 = sort(xm[xm < orio])
    t2 = sort(append(xm2,t2,after = length(xm2)))
    # if(length(t2) %% 3 == 0){
    #   par(mfrow = c(length(t2)%/%3,3))
    # }else{
    #   par(mfrow = c(length(t2)%/%3 + 1,3))
    # }
    d=data.frame(x=numeric(),y=numeric(),fill=character(),cluster=numeric())
    
    for(i in 1:length(t2)){
      ar = str_which(names(perlist),as.character(t2[i]))[1]
      output <- matrix(unlist(perlist[ar]), ncol = num_of_topics, byrow = FALSE)
      output=output[-nrow(output),]
      ## ggplot
      for (topic_list in 1:num_of_topics){  ## num_of_topics not found -- compute it
        ids=which(output[,topic_list]!=0)
        d=rbind(d,cbind(x=topic_list,y=as.numeric(output[ids,topic_list]),fill=let[ids],cluster=paste0("Cluster: ",t2[i],"-Elements: ",Clus$seqnum[t2[i]+1])))
      }
      d$fill=as.character(d$fill)
    }
    levels(d$fill)=let
    
    num_of_bins <- length(let)
    
    d$y=as.numeric(as.character(d$y))
    colors=topo.colors(num_of_bins)
    colors[num_of_bins]="cornsilk"
    colors[7]="yellow3"
    colors[1]="purple"
    colors[5]="seagreen"
    colors[6]="green"

    print({
    ggplot(data=d, aes(x=x, y=y, fill=fill)) +
      geom_bar(stat="identity") +
      #scale_fill_manual(values=heat.colors(num_of_bins))+
      scale_fill_manual(values= colors)+
      theme_bw() +
      #theme_minimal() + 
      facet_wrap(~cluster, scales='free', ncol = 3) +
      theme(legend.position="bottom") + theme(panel.spacing = unit(2, "lines"))
    })
    #ggsave(paste0(slash_for_topics,data_topic,"/","gg_bar_plots_level_",lev,".jpg"), width=5900, height=3900, res=100)
    
    dev.off() 
  }
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
}

Bar_plots_general <- function(g,jpeg_label,output_folder,slash_for_topics,data_topic,algo,num_of_topics){
  #if(flagtic == TRUE) tic(sprintf("Bar_plots_general -- Num of clusters: %d ", num_of_cl))
  jpeg(paste0(output_folder,slash_for_topics,data_topic,"/","gg_bar_plots","_algo_",algo,jpeg_label,".jpg"), width=3900, height=6000, res=100)
  
  d=data.frame(x=numeric(),y=numeric(),fill=character(),cluster=numeric())
  num_of_cl <- length(unique(g))
  for(i in unique(g)){
    a=data_discr[,] %>% filter(g==i)
    
    #metrics for evaluation 
    output <- create_matrix_for_evaluation(a,groups,let,sim)$permat
    output=output[-nrow(output),]
    ## ggplot
    for (topic_list in 1:num_of_topics){
      ids=which(output[,topic_list]!=0)
      d=rbind(d,cbind(x=topic_list,y=as.numeric(output[ids,topic_list]),fill=let[ids],cluster=paste0("Cluster: ",i,"-Elements: ",nrow(a))))
    }
    d$fill=as.character(d$fill)
  }
  levels(d$fill)=let
  d$y=as.numeric(as.character(d$y))
  colors=topo.colors(num_of_bins)
  colors[num_of_bins]="cornsilk"
  colors[7]="yellow3"
  colors[1]="purple"
  colors[5]="seagreen"
  colors[6]="green"
  
  print({
    ggplot(data=d, aes(x=x, y=y, fill=fill)) +
      geom_bar(stat="identity") +
      #scale_fill_manual(values=heat.colors(num_of_bins))+
      scale_fill_manual(values= colors)+
      theme_bw() +
      #theme_minimal() + 
      facet_wrap(~cluster, scales='free', ncol = 3) +
      theme(legend.position="bottom") + theme(panel.spacing = unit(2, "lines"))
  })
  #ggsave(paste0(slash_for_topics,data_topic,"/","gg_bar_plots_level_",lev,".jpg"), width=5900, height=3900, res=100)
  
  dev.off() 
  #if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
}

Bar_plots_general_without_output_folder <- function(g,jpeg_label){
  #if(flagtic == TRUE) tic(sprintf("Bar_plots_general -- Num of clusters: %d ", num_of_cl))
  jpeg(paste0(jpeg_label,".jpg"), width=3900, height=6000, res=100)
  
  d=data.frame(x=numeric(),y=numeric(),fill=character(),cluster=numeric())
  num_of_cl <- length(unique(g))
  for(i in unique(g)){
    a=data_discr[,] %>% filter(g==i)
    
    #metrics for evaluation 
    output <- create_matrix_for_evaluation(a,groups,let,sim)$permat
    output=output[-nrow(output),]
    ## ggplot
    for (topic_list in 1:num_of_topics){
      ids=which(output[,topic_list]!=0)
      d=rbind(d,cbind(x=topic_list,y=as.numeric(output[ids,topic_list]),fill=let[ids],cluster=paste0("Cluster: ",i,"-Elements: ",nrow(a))))
    }
    d$fill=as.character(d$fill)
  }
  levels(d$fill)=let
  d$y=as.numeric(as.character(d$y))
  colors=topo.colors(num_of_bins)
  colors[num_of_bins]="cornsilk"
  colors[7]="yellow3"
  colors[1]="purple"
  colors[5]="seagreen"
  colors[6]="green"
  
  print({
    ggplot(data=d, aes(x=x, y=y, fill=fill)) +
      geom_bar(stat="identity") +
      #scale_fill_manual(values=heat.colors(num_of_bins))+
      scale_fill_manual(values= colors)+
      theme_bw() +
      #theme_minimal() + 
      facet_wrap(~cluster, scales='free', ncol = 3) +
      theme(legend.position="bottom") + theme(panel.spacing = unit(2, "lines"))
  })
  #ggsave(paste0(slash_for_topics,data_topic,"/","gg_bar_plots_level_",lev,".jpg"), width=5900, height=3900, res=100)
  
  dev.off() 
  #if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
}

# A function to plot the barplot of a cluster
BarCl <- function(cl,perlist,let,flagtic,logFile){
  if(flagtic == TRUE) tic(sprintf("BarCl -- Cluster :%d ", cl))
  if(is.na(clep[cl])){
    print("Den yparxei")
  }else{
    #jpeg(paste0(output_folder,slash_for_topics,data_topic,"/","gg_bar_plot_cluster_",cl,".jpg"), width=5900, height=3900, res=100)
    #par(mar=c(3,3,4,4),xpd=TRUE)
    output <- matrix(unlist(perlist[cl+1]), ncol = str_length(udata$AA.JUNCTION[1]), byrow = FALSE)
    #barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d',cl)) 
    #legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)
    
    output=output[-nrow(output),]
    d=data.frame(x=numeric(),y=numeric(),fill=character())
    
    for (topic_list in 1:num_of_topics){
      ids=which(output[,topic_list]!=0)
      d=rbind(d,cbind(x=topic_list,y=as.numeric(output[ids,topic_list]),fill=let[ids]))
    }
    d$fill=as.character(d$fill)
    levels(d$fill)=let
    d$y=as.numeric(as.character(d$y))
    
    print({
    ggplot(data=d, aes(x=x, y=y, fill=fill)) +
      geom_bar(stat="identity") +
      #scale_fill_manual(values=heat.colors(num_of_bins))+
      scale_fill_manual(values= topo.colors(num_of_bins))+
      theme_bw() +
      #theme_minimal() + 
      theme(legend.position="bottom")
    })
    #dev.off() 
  }
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
}

# A function for showing  Sequences and Id's for specific level
AminoLev <- function(level,df,flagtic,logFile){
  if(flagtic == TRUE) tic(sprintf("AminoLev -- Level: %d ", level))
  if( is.element(lev,clep) == FALSE){
    print("Den yparxei")
  }else{
    sum(Clus[Clus$level == level,]$seqnum) # akoloy8ies sto level
    x4 = data_frame("Sequence.ID" = character(0),"AA.JUNCTION" = character(0))
    gg = na.omit(unique(df[,which(names(df) == sprintf("level.%d", level))]))
    for(i in 1:length(gg)){
      clust = gg[i]
      x3 = data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", clep[clust]))] == clust, ]$Sequence.ID), na.omit(df[df[which(names(df) == sprintf("level.%d", clep[clust]))] == clust, ]$AA.JUNCTION))
      names(x3) = c("Sequence.ID","AA.JUNCTION")
      se = x3$Sequence.ID
      aa = as.character(x3$AA.JUNCTION)
      aa = as.list(aa)
      names(aa) = se
      x4 = rbind(x4,x3)
    }
    se2 = x4$Sequence.ID
    aa2 = as.character(x4$AA.JUNCTION)
    aa2 = as.list(aa2)
    names(aa2) = se2
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    x4 # Print in console 
  }
}
#write.fasta(aa2, names= FALSE ,file = sprintf("Level.%d", level)) # Create a fasta file

# A function for showing  Sequences and Id's for specific cluster
AminoCl <- function(clust,clep,df,flagtic,logFile){
  #if(flagtic == TRUE) tic("AminoCl")
  #if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  if(clust == 0){
    x3 = data.frame("Sequence.ID" = df$Sequence.ID, "AA.JUNCTION" = df$AA.JUNCTION)
    return(x3)
  }else{
    if(!(is.na(clep[clust]))){
      x3 = data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", clep[clust]))] == clust, ]$Sequence.ID), na.omit(df[df[which(names(df) == sprintf("level.%d", clep[clust]))] == clust, ]$AA.JUNCTION))
      names(x3) = c("Sequence.ID","AA.JUNCTION")
      se = x3$Sequence.ID
      aa = as.character(x3$AA.JUNCTION)
      aa = as.list(aa)
      names(aa) = se 
      return(x3) # Print in console 
    }
  }
}

# A function which visualize the sequences of a data frame using absolute common letters (i.e. "A _ _ _ _ K R _ _ _ Q _ Y Y Y _ _ _ T _")
Opt <- function(df,flagtic,logFile){
  if(is.list(df)){
    #print("Den yparxei")
    #if(flagtic == TRUE) tic("Opt")
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      xar[i] = "_"
      if (f[1] == nrow(df)){
        xar[i] = names(f[1])
      }
    }
    # if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    paste(xar,collapse = ' ')
  }
}
#Opt(AminoCl(1))

# A function which visualize the sequences and percentages of a data frame using percentage common letters (i.e. "A _ _ _ _ K R _ _ _ Q _ Y Y Y _ _ _ T _")
OptNew <- function(df,logFile){
  if(is.list(df)){
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      if (length(which(f == max(f))) > 1){
        tempopt = "-"
        tempopt = paste(tempopt,names(which(f == max(f)))[1],sep = "")
        #tempopt =  names(which(f == max(f)))[1]
        for(j in 2:(length(which(f == max(f))))){
          tempopt = paste(tempopt, names(which(f == max(f)))[j],sep = "|")
        }
        tempopt = paste(tempopt,"-",sep = "")
        xar[i] = tempopt
      }else{
        xar[i] = names(which(f == max(f)))
      }
    }
    paste(xar,collapse = ' ')
  }
  
}
#OptNew(AminoCl(1))

# A function which visualize the sequences of a data frame using absolute similarity groups (i.e. "A _ _ _ _ K Am _ Ba _ Ac _ Y Y Y _ _ _ T _")
Opt2 <- function(df,alt,sim,altsim,flagtic,logFile){
  if(is.list(df)){
    #if(flagtic == TRUE) tic("Opt2")
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      xar[i] = "_"
      if (f[1] == nrow(df)){
        xar[i] = names(f[1])
      }else{
        y = TRUE
        d = str_which(sim[[i]],names(f[1]))
        for(j in 2:length(f)){
          y = y && str_detect(sim[[i]][d],names(f[j]))
        }
        if( y == TRUE){
          if(alt == TRUE){
            xar[i] = names(altsim[[i]][d])
          }else{
            xar[i] = names(sim[[i]][d])
          }
        } 
      }
      
    }
    #if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    paste(xar,collapse = ' ') 
  }
}
#Opt2(AminoCl(5),TRUE)

# A function which visualize the sequences and percentages of a data frame using absolute similarity groups (i.e. "A _ _ _ _ K Am _ Ba _ Ac _ Y Y Y _ _ _ T _")
Opt2New <- function(df,sim,logFile){
  if(is.list(df)){
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      for(j in 1:length(f)){
        if( length(which(str_detect(names(f),names(sim[str_which(sim,names(f)[j])])))) == 0 ){
          names(f)[j] = names(sim[str_which(sim,names(f)[j])])
        }else if(which(str_detect(names(f),names(sim[str_which(sim,names(f)[j])]))) != j){
          f[which(str_detect(names(f),names(sim[str_which(sim,names(f)[j])])))] = f[which(str_detect(names(f),names(sim[str_which(sim,names(f)[j])])))] + f[j]
          names(f)[j] = "dip"
        }
      }
      if(length(which(names(f) == "dip")) != 0 ){
        f = f[-which(names(f) == "dip")] 
      }
      if (length(which(f == max(f))) > 1){
        tempopt = "-"
        tempopt = paste(tempopt,names(which(f == max(f)))[1],sep = "")
        for(j in 2:(length(which(f == max(f))))){
          tempopt = paste(tempopt, names(which(f == max(f)))[j],sep = "|")
        }
        tempopt = paste(tempopt,"-",sep = "")
        xar[i] = tempopt
      }else{
        xar[i] = names(which(f == max(f)))
      }
    }
    paste(xar,collapse = ' ')
  }
}
#Opt2New(AminoCl(5))

# Returns the major #numOfTopics topics of input and the corresponding percentage of the non zero elements of each one of them.
findTopicsPerCl <- function(input,numOfTopics,let){
  if(is.list(input)){
    xar <- matrix(0,nrow=numOfTopics+1, ncol=str_length(input$AA.JUNCTION[1]))
    xar_num <- matrix(0,nrow=numOfTopics+1, ncol=str_length(input$AA.JUNCTION[1]))
    for (i in 1:str_length(input$AA.JUNCTION[1])) {
      f <- table(str_sub(input$AA.JUNCTION,i,i))
      a=sort(f,decreasing=T)
      if ((let[1] %in% names(a) & (a[let[1]]/nrow(input))>0.9) | (let[2] %in% names(a) & (a[let[2]]/nrow(input))>0.9) | 
          (let[1] %in% names(a) & let[2] %in% names(a) & (a[let[1]]/nrow(input)+a[let[2]]/nrow(input))>0.9)){
        #if (a[let[1]]/nrow(input)>0.9){
        xar[1:numOfTopics,i]=names(a[1:numOfTopics])
        xar[numOfTopics+1,i]=paste0(let[1],"|",let[2])
        xar_num[1:numOfTopics,i]=a[1:numOfTopics]
        xar_num[numOfTopics+1,i]=0
        if (let[1] %in% names(a)){
          xar_num[numOfTopics+1,i]=xar_num[numOfTopics+1,i]+a[let[1]]/nrow(input)}
        if (let[2] %in% names(a)){
          xar_num[numOfTopics+1,i]=xar_num[numOfTopics+1,i]+a[let[2]]/nrow(input)}
        #}
      }else{
        xar[numOfTopics+1,i]=paste0(let[1],"|",let[2])
        xar_num[numOfTopics+1,i]=0
        if (let[1] %in% names(a))
          xar_num[numOfTopics+1,i]=xar_num[numOfTopics+1,i]+a[let[1]]/nrow(input)
        if (let[2] %in% names(a))
          xar_num[numOfTopics+1,i]=xar_num[numOfTopics+1,i]+a[let[2]]/nrow(input)
        a=a[names(a)!=let[1]]
        a=a[names(a)!=let[2]]
        xar[1:numOfTopics,i]=names(a[1:numOfTopics])
        xar_num[1:numOfTopics,i]=a[1:numOfTopics]
      }
    }
    topics=1-xar_num[nrow(xar_num),][which(xar_num[nrow(xar_num),]!=1)] 
    names(topics)=which(xar_num[nrow(xar_num),]!=1)
    t=sort(topics,decreasing = T)
    t=t[1:numOfTopics]
  }
  return(t)
}

computeSimilarityCl<-function(per_matrix,letter_sim,let){
  new_freq_matrix=as.data.frame(matrix(0,nrow = nrow(per_matrix),ncol = ncol(per_matrix)))
  for (i in 1:ncol(per_matrix)){
    zero=which(per_matrix[,i]==0 & names(per_matrix[,i])!="Entropy")
    letter_freq_matrix=per_matrix[1:(nrow(per_matrix)-1),i]*letter_sim
    letter_freq_matrix[,zero]=0
    new_freq_matrix[1:(nrow(per_matrix)-1),i]=rowSums(letter_freq_matrix)/(length(let)-length(zero))
    new_freq_matrix[(nrow(per_matrix)),i]=entropy(new_freq_matrix[1:(nrow(per_matrix)-1),i], base=exp(1))
  }
  colnames(new_freq_matrix)=colnames(per_matrix)
  row.names(new_freq_matrix)=row.names(per_matrix)
  
  return(new_freq_matrix)
}

# A function plot a diagramm with Identity or Similarity of clusters in every level
Id <- function(ff,cho,flagtic,logFile){
  if(flagtic == TRUE) tic(sprintf("Id -- Type: $s ", cho))
  if(cho == "Identity"){
    par(xpd=TRUE)
    pp <- data.frame(x = ff$level,y = ff$sumper)
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggplot(na.omit(pp),aes(x = x,y = y)) + stat_sum()
  }else{
    par(xpd=TRUE)
    pp <- data.frame(x = ff$level,y = ff$sumper2)
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggplot(na.omit(pp),aes(x = x,y = y)) + stat_sum()
  }
}

# A function plots a collapseble tree with id, sequences, identity and similarity of every cluster - node 
cc <- function(df,Clus,flagtic,logFile){
  if(flagtic == TRUE) tic()
  df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
  ii = Clus$seqnum
  ii = as.character(ii)
  
  trid1 = df[str_which(names(df), "level.")] # Data frame with identities percentage
  kk1 = trid1
  trsim = df[str_which(names(df), "level.")] # Data frame with similarities percentage
  iii = df[str_which(names(df), "level.")]
  for(i in 1:length(trid1)){
    tempg = trid1[,i]+1
    trid1[i] = round(Clus$Identity[tempg],digits = 2)
    temph = trsim[,i]+1
    trsim[i] = round(Clus$Similarity[temph], digits = 2)
    tempi = iii[,i] + 1
    iii[i] = ii[tempi]
  }
  
  for(i in 2:length(str_which(names(df), "level."))){
    trid1[,i] = paste(kk1[,i],iii[,i],trid1[,i], trsim[,i],sep = " ")
    for(j in 1:nrow(df)){
      if(str_detect(trid1[j,i],"NA") == TRUE){
        trid1[j,i] = NA
      }
    }
  }
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0("Collapsible Tree","\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
  collapsibleTree(
    trid1,
    hierarchy = colnames(df[str_which(names(df), "level.")]),
    fill = c("jj",na.omit(ii)),
    width = 1820,
    height = 775,
    collapsed = FALSE
  )
}

# A function show the value of identity or similarity for level clusters and leaves until this level
idenlev <- function(lev,clep,Clus,xN,cho,flagtic,logFile){
  if(flagtic == TRUE) tic(sprintf("idenlev -- Level: %d ", lev))
  if(cho == "Identity"){
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    t1 = which(clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    sum(na.omit(Clus[t1+1,]$Identity)) /  length(na.omit(Clus[t1+1,]$Identity))
  }else{
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    t1 = which(clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    sum(na.omit(Clus[t1+1,]$Similarity)) /  length(na.omit(Clus[t1+1,]$Similarity))
  }
}

# A function show the value of identity or similarity for a cluster
idencl <- function(cl,Clus,cho,flagtic,logFile){
  if(flagtic == TRUE) tic(sprintf("idencl -- Cluster :%d ", cl))
  if(is.element(cl,Clus$ClusterId)){
    if(cho == "Identity"){
      if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
      na.omit(Clus[Clus$ClusterId == cl,])$Identity
    }else{
      if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
      na.omit(Clus[Clus$ClusterId == cl,])$Similarity 
    }
  }else{
    "Den yparxei"
  } 
}

# A function create a table with average identity, identity standard deviation, average similarity and similarity standard deviation for every level 
EmPin <- function(Clus,clep,xN,dfsd,flagtic,logFile){
  if(flagtic == TRUE) tic("EmPin")
  print(max(na.omit(clep)))
  for (i in 0:(max(na.omit(clep)))) {
    if(i == 0){
      t0=1
       
      m <- matrix(ncol = length(evaluation_metrics), nrow = 1, 0)
      for (j in 3:length(evaluation_metrics)){
        calculation <- strsplit(evaluation_metrics[j],"-")[[1]][1]
        metric <- strsplit(evaluation_metrics[j],"-")[[1]][2]
        if (calculation=="Average"){
          m[j] <- mean(Clus[[metric]][i+1])
        }else{
          m[j] <- sd(Clus[[metric]][i+1])
        }
      }
      
    }else{
      t = which(clep == i)
      # epipleon
      xm = as.numeric(as.data.frame(xN$leaves))
      orio = min(which(clep == i))
      xm2 = sort(xm[xm < orio])
      t = sort(append(xm2,t,after = length(xm2)))
      
      t0=length(t)
      
      m <- matrix(ncol = length(evaluation_metrics), nrow = 1, 0)
      for (j in 3:length(evaluation_metrics)){
        calculation <- strsplit(evaluation_metrics[j],"-")[[1]][1]
        metric <- strsplit(evaluation_metrics[j],"-")[[1]][2]
        if (calculation=="Average"){
          m[j] <- mean(na.omit(Clus[[metric]][t+1]))
        }else{
          m[j] <- sd(na.omit(Clus[[metric]][t+1]))
        }
      }
    }
    rows = sprintf("level.%d", i)
    dfsd[i+1,1] = c(rows)
    dfsd[i+1,2:ncol(dfsd)] = c(t0,as.vector(m)[3:ncol(m)])
  }
  colnames(dfsd)[1]="Level"
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  dfsd
}

normalize <- function(x) {
  return ((x - min(x[is.na(x) == FALSE])) / (max(x[is.na(x) == FALSE]) - min(x[is.na(x) == FALSE])))
}


###### Graph 
computeDistances <- function(xN,lev,df,Clus,use_only_leafs,sim,altsim,clep,num_of_topics){
  if(flagtic == TRUE) tic(sprintf("computeDistances -- level: %d", lev))
  
  print(paste0("compute distances level ", lev))
  if (use_only_leafs){
    t1 = which(clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(clep == lev))
    xm2 = sort(xm[xm < orio])
    useful_nodes = sort(append(xm2,t1,after = length(xm2)))
    bo=length(useful_nodes)
  }else{
    bo = length(which(!is.na(clep)))
    useful_nodes=which(!is.na(clep))
  }
  
  
  ffg = as.data.frame(matrix(0,nrow = bo,ncol = (bo))) # the distance between 2 Nodes
  ffg2 = as.data.frame(matrix(0,nrow = bo,ncol = (bo)))
  ffg2sim = as.data.frame(matrix(0,nrow = bo,ncol = (bo)))
  ffg3 = as.data.frame(matrix(NA,nrow = bo,ncol = (bo))) # the final distance combining ffg and ffg2
  
  print("start with sapply")
  print(bo)
  s = sapply(1:(bo), function(i){
    tempN1 = FindNode(xN,(sprintf("%d", useful_nodes[i])))
    tempN1 = tempN1$path
    temp2N1 = Opt(AminoCl(useful_nodes[i],clep,df,flagtic,logFile),flagtic,logFile)
    temp3N1 = Opt2(AminoCl(useful_nodes[i],clep,df,flagtic,logFile),TRUE,sim,altsim,flagtic,logFile)
    sapply(i:(bo),function(j) {
      tempN2 = FindNode(xN,(sprintf("%d", useful_nodes[j])))
      tempN2 = tempN2$path
      tem = which(tempN1 == tempN2)
      d1 = length(tempN1) - tem[length(tem)]
      d2 = length(tempN2) - tem[length(tem)]
      temp2N2 = Opt(AminoCl(useful_nodes[j],clep,df,flagtic,logFile),flagtic,logFile)
      temp2N2 = str_replace_all(temp2N2,"_"," ")
      temp3N2 = Opt2(AminoCl(useful_nodes[j],clep,df,flagtic,logFile),TRUE,sim,altsim,flagtic,logFile)
      temp3N2 = str_replace_all(temp3N2,"_"," ")
      set(ffg, i, j, as.integer(d1 + d2)) 
      set(ffg2, i, j, stringdist(temp2N1,temp2N2,method = "lv" )) 
      set(ffg2sim, i, j, stringdist(temp3N1,temp3N2,method = "lv" )) 
    })
  })
  
  #i = 1
  # anw trigwnikos
  
  #find max
  ma = max(ffg[is.na(ffg) == FALSE])
  ffg[ffg == 0] = ma
  ffg2[ffg2 == 0] = num_of_topics
  ffg2sim[ffg2sim == 0] = num_of_topics
  ffg[lower.tri(ffg)] = NA
  ffg2[lower.tri(ffg2)] = NA
  ffg2sim[lower.tri(ffg2sim)] = NA
  ffg = ffg / ma * 100
  ffg2 = num_of_topics - ffg2
  ffg2 = ffg2 / num_of_topics * 100
  ffg2sim = num_of_topics - ffg2sim
  ffg2sim = ffg2sim / num_of_topics * 100
  
  diag(ffg) = 0
  diag(ffg2) = 0
  diag(ffg2sim) = 0
  diag(ffg3) = 0
  ffg3 = (ffg2 + ffg2sim + ffg) / 3
  
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  listnet = list("ffg" = ffg, "ffg2" = ffg2, "ffg2sim" = ffg2sim, "ffg3" = ffg3)

  return(listnet)
}

Netw <- function(xN,clep,lev,thr,thrt,netyp,df,Clus,use_only_leafs, ffg, ffg2, ffg3, ffg2sim, net_sil){
  #if(flagtic == TRUE) tic(sprintf("Network -- level: %d, Network type: %s, Threshold type: %s, Threshold value: %d, Using silhouette: %d ", lev,netyp,thrt,thr,net_sil))
  if (use_only_leafs){
    t1 = which(clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(clep == lev))
    xm2 = sort(xm[xm < orio])
    useful_nodes = sort(append(xm2,t1,after = length(xm2)))
    bo=length(useful_nodes)
  }else{
    bo = length(which(!is.na(clep)))
    useful_nodes=which(!is.na(clep))
  }
  
  if(thrt == "Distance"){
    thrtyp = ffg
  }else if (thrt == "StrIdentSimilarity"){
    thrtyp = ffg2
  }else if (thrt == "StrGroupSimilarity"){
    thrtyp = ffg2sim
  }else{
    thrtyp = ffg3
  }
  
  tempor = ffg
  tempor2 = ffg2
  tempor2sim = ffg2sim
  tempor3 = ffg3
  tempor[thrtyp > thr] = 0
  tempor2[thrtyp > thr] = 0
  tempor2sim[thrtyp > thr] = 0
  tempor3[thrtyp> thr] = 0 #when the distance between two nodes is too big then make it equal to 0 in order to delete the edge latter

  if(thrt == "Distance"){
    thrtyp2 = tempor
  }else if (thrt == "StrIdentSimilarity"){
    thrtyp2 = tempor2
  }else if (thrt == "StrGroupSimilarity"){
    thrtyp2 = tempor2sim
  }else{
    thrtyp2 = tempor3
  }
  
  if(net_sil == TRUE & !use_only_leafs){
    newffg = normalize(ffg)
    newffg2 = normalize(ffg2)
    newffg2sim = normalize(ffg2sim)
    newffg3 = newffg * 0.5 + newffg2 * 0.25 + newffg2sim * 0.25  
    nnnn = normalize(newffg3)
    
    jhj = silhouette(c(0,Clus$level[useful_nodes+1]),t(nnnn))
    matches <- regmatches(unique(x$pathString), gregexpr("[[:digit:]]+", unique(x$pathString)))
    un=c(0,useful_nodes)
    tttsyn = lapply(1:length(matches),function(i){
      m=c()
      for (k in 1:length(as.numeric(unlist(matches[i])))){
        print(which(un==as.numeric(unlist(matches[i]))[k]))
        m=c(m,which(un==as.numeric(unlist(matches[i]))[k]))
      }
      ttt = jhj[m,3]
      sapply(length(ttt):1,function(j){
        if(j > 1){
          if(ttt[j] >= ttt[j-1]){
            set(tempor3, j-1, j, 0)
            set(ffg3, j-1, j, 0)
            set(thrtyp, j-1, j, 0)
            set(thrtyp2, j-1, j, 0)
          }
        }
      })
    })
  }
  
  jjj = matrix(1,nrow = (bo),ncol = ((bo)))
  #jjj = matrix(1,nrow = length(bb),ncol = length(bb))
  net0 <- graph_from_adjacency_matrix(jjj,mode = "upper")
  net1 <- graph_from_adjacency_matrix(jjj,mode = "upper")
  
  bb = vector(length = (length(useful_nodes)))
  #bb[1]=0
  bb[1:length(bb)] = clep[useful_nodes]
  #XRWMA
  # Generate colors based on media type:
  colrs <- c("#1E90FF", "#BA55D3", "#0000FF", "#557fd2", "#54d17e", "#8aad62", "#C6E2FF", "#e5e234", "#FFD700", "#00EE00", "#C1FFC1", "#ea8509", "#54FF9F", "#FF0000", "#ed3b1c", "#ed1c7a", "#0c0c0c", "#b8d8af", "#ED9121","#45f713")
  colrs = colrs[1:(max(na.omit(bb)))]
  V(net0)$color <- colrs[bb]
  V(net1)$color <- colrs[bb]
  
  bb = vector(length = (max(bo)))
  #bb[1]=0
  bb[1:length(bb)] = useful_nodes
  
  # Compute node degrees (#links) and use that to set node size:
  deg <- (Clus$seqnum[bb+1] / (max(Clus$seqnum[bb+1]))) * 20
  #deg[1]=max(deg[2:length(deg)])
  V(net0)$size <- deg
  V(net1)$size <- deg
  
  # The labels are currently node IDs.
  # Setting them to NA will render no labels:
  V(net0)$label <- Clus$ClusterId[bb+1] # or V(net0)$label <- useful_nodes
  V(net1)$label <- Clus$ClusterId[bb+1]
  
  # Set edge width based on weight:
  hhh = na.omit(as.vector(t(ffg3))) #vector of number of elements equal to bo x bo
  hhh1 = na.omit(as.vector(t(tempor3))) #vector of number of elements equal to bo x bo
  E(net0)$width <- hhh 
  E(net1)$width <- hhh1 
  E(net0)$weight <- hhh 
  E(net1)$weight <- hhh1 
  
  #change arrow size and edge color:
  E(net0)$arrow.size <- .2
  E(net0)$edge.color <- "gray80"
  E(net1)$arrow.size <- .2
  E(net1)$edge.color <- "gray80"
  pal1 <- rainbow(6, alpha=1) 
  
  net0.copy <- igraph::delete.edges(net0, which(E(net0)$width == 0))
  #net0.copy <- igraph::delete.edges(net0, is.na(V(net0)$size))
  net0.copy <- igraph::delete.vertices(net0.copy, is.na(V(net0)$size))
  net1.copy <- igraph::delete.edges(net1, which(E(net1)$width == 0))
  
  colors=c("red","blue","green","pink","grey") #red:strong relationship
  #colors=c("grey","pink","green","red","blue") 
  w=c(2.5,2,1.5,1,0.5)/5
  
  step1=(max(tempor3, na.rm = T)-min(tempor3, na.rm = T))/5
  
  for (i in 1:5){
    E(net1.copy)[weight > (i*step1-step1) & weight<i*step1]$color <-colors[i]
    E(net1.copy)[weight > (i*step1-step1) & weight<i*step1]$width <-w[i]
  }
  E(net1.copy)[weight == 0]$color <-colors[5] #0 represents to 100% dissimilarity
  E(net1.copy)[weight== 0]$width <-w[5]
  
  step=(max(ffg3, na.rm = T)-min(ffg3, na.rm = T))/5
  E(net0.copy)$color <-"grey"
  
  for (i in 1:5){
    E(net0.copy)[weight > (i*step1-step1) & weight<i*step1]$color <-colors[i]
    E(net0.copy)[weight > (i*step-step) & weight<i*step]$width <-w[i]
  }
  E(net0.copy)[weight == 0]$color <-colors[5]
  E(net0.copy)[weight== 0]$width <-w[5]

  # %/%: integer division
  
  plot(net0.copy, edge.curved=.1, vertex.label.color = "black", edge.width=1)  #plot the network graph
  legend("topleft", inset=c(0.1,0.2), colnames(df[str_which(names(df), "level.")])[1:(lev+1)], pch=21,
         col="#777777", pt.bg=unique(V(net0)$color), pt.cex=2, cex=.8, bty="n", ncol=1)
  legend("topright", inset=c(0.1,0.2), c("4-5","3-4","2-3","1-3","0-1"), pch=21,
         col="#777777", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
  #matrix.heatmap(thrtyp)
  
  #relationship strength=5 means that the nodes are very similar to each other 
  
  plot(net1.copy, edge.curved=.1, vertex.label.color = "black",edge.width=2)  #plot the network graph
  legend("topleft", inset=c(0.1,0.2), paste(unique(c(0,clep[useful_nodes]))), pch=21,
         col="#777777", pt.bg=unique(V(net1)$color), pt.cex=2, cex=.8, bty="n", ncol=1)
  legend("topright", inset=c(0.1,0.2), c("4-5","3-4","2-3","1-3","0-1"), pch=21,
         col="#777777", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
  
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  listnet = list("net0.copy" = net0.copy, "net1.copy" = net1.copy, "bb" = bb, "colrs" = colrs, "tempor" = tempor, "tempor2" = tempor2, "tempor2sim" = tempor2sim, "tempor3" = tempor3,"thrtyp" = thrtyp, "thrtyp2" = thrtyp2)
  return(listnet)
}


##### Evaluation 
create_matrix_for_evaluation <- function(data,groups,let,sim){
  data_new=data
  
  #Discretize data
  for (i in 1:(length(groups)-1)){
    ids=which(data[,2:ncol(data)]>=groups[i] & data[,2:ncol(data)]<groups[i+1],arr.ind=TRUE)
    for (j in unique(ids[,2]))
      data_new[ids[which(ids[,2]==j),1],j+1]=let[i]
    if (i==length(groups)-1){
      ids=which(data[,2:ncol(data)]==groups[i+1],arr.ind=TRUE)
      for (j in unique(ids[,2]))
        data_new[ids[which(ids[,2]==j),1],j+1]= let[i] #let[i+1]
    }
  }
  
  data_new$x <- apply( data_new[ , 2:ncol(data_new) ] , 1 , paste , collapse = "" )
  udata <- data_new[c(1,ncol(data_new))]
  colnames(udata)=c("Sequence.ID","AA.JUNCTION")
  
  #Create frequency matrices 
  mymat = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
  permat = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  simmat = matrix(0,nrow = max_group_length,ncol =str_length(udata$AA.JUNCTION[1]))
  persim = matrix(0,nrow = max_group_length+1,ncol =str_length(udata$AA.JUNCTION[1]))
  
  rownames(mymat) = let
  rownames(permat) = c(let,"Entropy")
  rownames(simmat) = c(paste0("group_",1:max_group_length))
  rownames(persim) = c(paste0("group_",1:max_group_length),"Entropy")
  
  trimudata = strsplit(udata$AA.JUNCTION,"")
  align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
  
  for(i in 1:str_length(udata$AA.JUNCTION[1])){
    temptab = plyr::count(align[i],vars = colnames(align)[i])
    names(temptab)[1] = "X1"
    match_elements=match(names(mymat[,i]),as.vector(unlist(temptab[1])))
    ids=as.vector(unlist(temptab[1]))[na.omit(match_elements)]
    mymat[ids,i] = as.vector(unlist(temptab[2]))[na.omit(match_elements)]
    #mymat[which(is.na(match(names(mymat[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
    permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
    temptab1 = temptab
    temptab1[1] = names(sim[[i]])[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim[[i]],x)))]
    temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
  }
  permat[1:length(let),] = (mymat / length(udata$AA.JUNCTION)) * 100
  persim[1:max_group_length,] = (simmat / length(udata$AA.JUNCTION)) * 100
  
  identity=length(which(permat==100))/str_length(udata$AA.JUNCTION[1])
  similarity=length(which(persim==100))/str_length(udata$AA.JUNCTION[1])
  entropy_id=mean(permat[nrow(permat),])
  entropy_sim=mean(persim[nrow(persim),])
  
  result = list("udata" = udata,"permat"= permat, "persim" = persim)
  return(result)
}

create_freq_matrix_for_merged_clusters <- function(udata,groups,let,sim,max_group_length){
  #Create frequency matrices 
  mymat = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
  permat = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  simmat = matrix(0,nrow = max_group_length,ncol =str_length(udata$AA.JUNCTION[1]))
  persim = matrix(0,nrow = max_group_length+1,ncol =str_length(udata$AA.JUNCTION[1]))
  
  rownames(mymat) = let
  rownames(permat) = c(let,"Entropy")
  rownames(simmat) = c(paste0("group_",1:max_group_length))
  rownames(persim) = c(paste0("group_",1:max_group_length),"Entropy")
  
  trimudata = strsplit(udata$AA.JUNCTION,"")
  align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
  
  for(i in 1:str_length(udata$AA.JUNCTION[1])){
    temptab = plyr::count(align[i],vars = colnames(align)[i])
    names(temptab)[1] = "X1"
    mymat[which(is.na(match(names(mymat[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
    permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
    temptab1 = temptab
    gg = as.vector(unlist(temptab1[1]))
    temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim[[i]])[which(str_detect(unlist(sim[[i]]),gg[x]))]))) #delete digits from the names of similarity groups
    temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
    simmat[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),i] = as.vector(unlist(temptab1[2]))
    persim[max_group_length+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
    
  }
  permat[1:length(let),] = (mymat / length(udata$AA.JUNCTION)) * 100
  persim[1:max_group_length,] = (simmat / length(udata$AA.JUNCTION)) * 100
  
  identity=length(which(permat==100))/str_length(udata$AA.JUNCTION[1])
  similarity=length(which(persim==100))/str_length(udata$AA.JUNCTION[1])
  entropy_id=mean(permat[nrow(permat),])
  entropy_sim=mean(persim[nrow(persim),])
  
  result = list("udata" = udata,"permat"= permat, "persim" = persim)
  return(result)
}


##### Topic Similarity 
compute_topic_similarity <- function(){
  ########### Compute Topic Similarity ########### 
  words_per_topic=20
  topic_words=read.csv(paste0(dataset,slash_for_topics,data_topic,"/LDA/model-final.twords"),sep=" ",header = F,stringsAsFactors = F)
  #Sys.setlocale(category = "LC_ALL", locale = "Greek")
  #topic_words=fix(topic_words)
  
  t=topic_words %>% group_by(V1) %>% summarise(n=n())
  t=t[order(-t$n),]
  t=t %>% filter(t$n>1 & V1!="Topic")
  
  topic_words_list=c()
  
  row=-20+1
  for (i in 1:num_of_topics){
    row=row+20+1
    topic_words_list=cbind(topic_words_list,topic_words[row:(row+20-1),1]) 
  }
  topic_words_list=as.data.frame(topic_words_list,stringsAsFactors=F)
  colnames(topic_words_list)=paste0("topic_",1:num_of_topics)
  
  topic_similarity = as.data.frame(matrix(0,nrow = num_of_topics,ncol = num_of_topics))
  
  for (i in 1:num_of_topics){
    for (j in 1:num_of_topics){
      #if (j>i){
      all_pr=data.frame(V=c(topic_words_list[,i],topic_words_list[,j]))
      all_pr=all_pr %>% group_by(V) %>% summarise(n=n()) %>% filter(n>1)
      topic_similarity[i,j]=nrow(all_pr)/20
      #}
    }
  }
}


##### Semantic Similarity
getGoogleCount <- function(searchTerms=NULL, language="en", ...){
  require(RCurl)
  entry    <- paste(searchTerms, collapse="+")
  siteHTML <- getForm("http://www.google.com/search",
                      hl=language, lr="", q=entry,
                      btnG="Search")
  
  write.table(siteHTML, file="tmp google.txt")  
  indicatorWord <- "resultStats"        
  posExtractStart <- gregexpr(indicatorWord, siteHTML,
                              fixed = TRUE)[[1]]
  stringExtract <- as.character(substring(siteHTML, first=posExtractStart[2]-30,
                                          last = posExtractStart[2] +50 ))
  count <- strsplit(stringExtract, 'resultStats')[[1]][2] 
  count <- strsplit(count, split='results')[[1]][1]
  count <- strsplit(count, split='>')[[1]][2]
  if(length(strsplit(count, split=" ")[[1]])==2){
    count <- strsplit(count, split=" ")[[1]][2] 
  }
  count <- as.numeric(gsub(",", "", count))
  return(count)
}

getDistance <- function(x,y){
  xy <- getGoogleCount(c(x, y)) 
  x  <- getGoogleCount(c(x))
  y  <- getGoogleCount(c(y))
  
  xy <- as.numeric(gsub(",", "", xy))
  x  <- as.numeric(gsub(",", "", x ))
  y  <- as.numeric(gsub(",", "", y ))
  M <- 859000000 
  dist <- (max(log(x), log(y)) - log(xy))/(log(M)-min(log(x), log(y)))   
  return(dist)
}

NGD <- function(words, language="en", print=FALSE,list=FALSE, ...){
  
  # check for arguments
  #if(!hasArg(words)) stop('NGD needs TWO strings like
  #                        c("word","word2") as word argument!')
  if(length(words)!=2) stop('word arguments has to be of
                            length two, e.g. c("word","word2")')
  
  # M: total number of web pages searched by google (2007)
  if(hasArg(M)) M <- list(...)$M else M <- 8058044651    
  
  x <- words[1]
  y <- words[2]
  
  # using getGoogleCount() function (see here)
  freq.x  <- getGoogleCount(x, language=language)
  freq.y  <- getGoogleCount(y, language=language)
  freq.xy <- getGoogleCount(c(x,y), language=language)
  
  # apply formula
  NGD = (max(log(freq.x), log(freq.y)) - log(freq.xy)) /
    (log(M) - min( log(freq.x), log(freq.y)) )
  
  # print results to console if requested
  if(print==TRUE){
    cat("\t", x,":", freq.x, "\n",
        "\t", y,":", freq.y, "\n",
        "\t", x,"+", y,":", freq.xy, "\n",
        "\t", "normalized google distance (NGD):",
        NGD, "\n", "\n")
  }
  
  
  # return list of results if requested (no default)
  # containing NGD and all counts. As default only one
  # the NGD is returned as numeric value
  
  results <- list(NGD=NGD,
                  x=c(x, freq.x),
                  y=c(y, freq.y),
                  xy=c(paste(x,"+",y), freq.xy)) 
  
  if(list==TRUE) return(results) else  return(NGD)
}
