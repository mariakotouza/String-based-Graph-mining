#library(webshot); #webshot::install_phantomjs() #in case phantomjs was not installed 
library("NMF")
library("NMI")
#library("mstknnclust")
library(rfUtilities)

source("functions.R")
source("global.R")
#source("SSL.R")

# for ssl
# source("https://bioconductor.org/biocLite.R")
# biocLite("graph")
# install.packages("NetPreProc")

data <- data.frame()
values <- list(flagEX=0,filterdf=data.frame(Columns=character(),Keys=character(),"I/E"=character()),ig=NULL,indexes=c(),forest=c(),sim=NULL,uni=NULL,uniframe=NULL)
clusterValues <- list()
mstValues<- list(edges=NULL,centroids=NULL,clusters=NULL,flag2=NULL,legenddf=NULL,legenddf2=NULL,shape=NULL)
centralValues <- list()
clustering_eval_matrix <- data.frame()

run_graph_tool <- function(dataset_name, datapath, sep, header, pipeline, group_by_clusterId, seqSelect, simSelect, 
                           componentSelect, excludeButton, includeButton, graph.type,
                           select1, slider1, text1, mstAlgoSelect, delete_edges_of_diff_classes, clusterMST, colormst, 
                           bordermst, fastCButton, centralSelect, clusterCentral, 
                           centralColor, clustering_algo_Select, clusterColor, 
                           plotHeat, heatSelect, percentage_labeled_data, iter) {
  
  if (save_tables_individually){
    #output folder 
    if(!file.exists(paste0(tmp_path,"/output"))){ 
      dir.create(paste0(tmp_path,"/output"))
    }
    if(!file.exists(paste0(tmp_path,"/output/", simSelect))){ 
      dir.create(paste0(tmp_path,"/output/", simSelect))
    }
    if(!file.exists(paste0(tmp_path,"/output/", simSelect, "/", dataset_name))){ 
      dir.create(paste0(tmp_path,"/output/",simSelect, "/", dataset_name))
    }
    output_folder <<- paste0(tmp_path,"/output/", simSelect, "/", dataset_name)
    
  }
  
  k_knn <- NA
  graph_thr <- NA
  pipeline <- strsplit(pipeline, ",")[[1]]
  
  clusterSelect1 <- clustering_algo_Select[1]
  clusterSelect2 <- clustering_algo_Select[2]
  
  read_data(datapath, sep="\t", header=T)
  calculate_distance(seqSelect, simSelect)
  all_ids <- 1:nrow(data)
  
  if (percentage_labeled_data == 0){
    data$V2 <<- NA 
  }else if (percentage_labeled_data < 1){
    # Take known labels equaly from each class
    known.label <- c()
    actual_classes <- unique(data$V1)
    for (k in 1:length(actual_classes)){
      s_ids <- which(data$V1==actual_classes[k])
      known.label <- c(known.label, sample(s_ids, round(percentage_labeled_data*length(s_ids))) )
    }
    
    #known.label <- sample(all_ids, percentage_labeled_data*length(all_ids))
    data$V2 <<- NA 
    data$V2[known.label] <<- data$V1[known.label]

    print(paste("Percentage of known labels", length(which(!is.na(data$V2[known.label])))/length(all_ids) ))
  }else{
    data$V2 <<- data$V1
  }
  
  print(unique(data$V2))
  
  if (1 %in% pipeline){
    filter_data(excludeButton, includeButton, select1, slider1, text1, ReverseButton)
  }
  
  if (length(graph.type)>0){
    for (type in graph.type){
      temp <- strsplit(type,":")[[1]]
      graph.type <- temp[1]
      if ((temp[1]=="knn") | (temp[1]=="mknn")){
        k_knn <- as.numeric(temp[2])
      }else{
        graph_thr <- as.numeric(temp[2])
      }
    }
  }
  
  ig<- list()
  for (type in graph.type){
    ig[[type]] <- graph_construction(clusterId=group_by_clusterId, seqSelect=seqSelect, graph.type ,k=k_knn, alpha=1, alpha1=-2,alpha2=1,
                           graph_thr, delete_edges_of_diff_classes=delete_edges_of_diff_classes)
  }
  
  if (3 %in% pipeline){
    select_component(componentSelect)
  }
  
  adj_table()
  
  if (2 %in% pipeline){
    visualize_network(graph_thr) 
  }
  
  download_graph() 

  #neigh_table() 
  
  #if (delete_edges_of_diff_classes){
    with_edges <- paste0("classes_", percentage_labeled_data)
  #}else{
    #with_edges <- "with_diffclass_edges"
  #}
  
  if ((graph.type=="knn") | (graph.type=="mknn")){
    thr <- k_knn
  }else{
    thr <-graph_thr  
  }
  
  #Find the percentage of data points that have been assigned to the correct cluster
  clustering_eval_matrix <- as.data.frame(matrix(0, nrow = length(clustering_algo_Select), ncol = 11))
  colnames(clustering_eval_matrix) <- c("Algorithm", "N_actual", "N", "Fmeasure", "Silhouette", "Purity", "Entropy", 
                                        "NMI","Modularity", "Conductance", "Time")
  
  #Clustering
  k <- 0
  if (4 %in% pipeline){
    for (cl in clustering_algo_Select){
      k <- k + 1
      t <- clustering(cl)
      if (t!=-1){
        ss <- silhouette(clusterValues$membership[[cl]], cl)
        clustering_eval_matrix$Silhouette[k] <- ss
        
        # The major class of each cluster
        clusterValues$membership[[cl]]
        data$V1
        
        actual_classes <- data$V1
        labels <- unique(actual_classes)
        g2 <- clusterValues$membership[[cl]]
        # Create a matrix with all the clustering methods and the metrics
        
        # F-measure
        Fscore <- 0
        F2 <- c()
        num_of_objects_class_r <- c()
        pres_all <- c()
        
        for (r in 1:length(unique(actual_classes))){
          # find the documents that belong to class r and save them into dataframe content 
          nr <- length(which(actual_classes==labels[r]))
          FLr <- c()
          pres_all_temp <- c()
          recall_all_temp <- c()
          #for (dif_levels in 1:length(unique(actual_classes))){
          #g2 <- cutree(hClusters_diana, k=dif_levels)
          for (i in unique(g2)){
            #find the documents that belong to cluster i
            content <- which(g2==i)
            ni <- length(content)
            nri <- length(which(actual_classes[content]==labels[r]))
            precision <- nri/ni
            recall <- nri/nr
            if ((precision+recall)>0){
              F1 <- 2*precision*recall/(precision+recall)
            }else{
              F1 <- 0
            }
            FLr <- c(FLr, F1)
            pres_all_temp <- c(pres_all_temp,precision)
            recall_all_temp <- c(recall_all_temp,recall)
          }
          #}
          
          #for each class find the maximum Fscore from the tree levels 
          F2 <- c(F2, max(FLr, na.rm = T))
          num_of_objects_class_r <- c(num_of_objects_class_r, nri)
          pres_all <- c(pres_all, max(pres_all_temp))
          Fscore <- Fscore + max(FLr, na.rm = T) * nr/nrow(data)
        }
        F2 <- sum(F2*num_of_objects_class_r)/sum(num_of_objects_class_r)
        #print(length(as.factor(data$V1)))
        #print(i)
        #print(((clusterValues$membership[[i]])))
        clustering_eval_matrix$Fmeasure[k] <- F2
        
        clustering_eval_matrix$Algorithm[k] <- cl
        clustering_eval_matrix$Time[k] <- t
        
        # Entropy 
        clustering_eval_matrix$Entropy[k] <- NMF::entropy(as.factor(data$V1), as.factor(clusterValues$membership[[cl]]))
        
        # Purity 
        clustering_eval_matrix$Purity[k] <- NMF::purity(as.factor(data$V1), as.factor(clusterValues$membership[[cl]]))
        clustering_eval_matrix$N_actual[k] <- length(unique(actual_classes))
        clustering_eval_matrix$N[k] <- length(unique(clusterValues$membership[[cl]]))
        
        # Nmimeasure
        tryCatch({
          clustering_eval_matrix$NMI[k] <- aricode::NMI(data$V1, clusterValues$membership[[cl]])
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        # Graph metrics
        clustering_eval_matrix$Modularity[k] <- modularity(ig[[type]],clusterValues$membership[[cl]])
        clustering_eval_matrix$Conductance[k] <- conductance(ig[[type]],clusterValues$membership[[cl]])$conductance
        
        if (5 %in% pipeline){
          visualize_clustered_graph(clusterColor, cl)
        }
      }
      
      
    }
    if(!file.exists(paste0(output_folder,"/clustering"))){ 
      dir.create(paste0(output_folder,"/clustering"))
    }
    write.table(clustering_eval_matrix, file = paste0(output_folder,"/clustering", "/", "clustering_", 
                                                      with_edges, "_", graph.type, "_", thr,"_",iter, ".txt"), row.names = F)
    print(clustering_eval_matrix)
  }
  
  if (6 %in% pipeline){
    #clustering_heatmap(plotHeat)
    #silhouette() 
    clustering_common(clusterSelect1, clusterSelect2)
    clustering_confusionMatrix1(clusterSelect1, clusterSelect2)
    clustering_confusionMatrix2(clusterSelect1, clusterSelect2)
    clustering_metrics() 
    #clustering_interheat(heatSelect)
  }
  
  #MST
  if (7 %in% pipeline){
    create_mst(mstAlgoSelect) 
    adj2_table() 
    ig_mst <- mst_construction(colormst, bordermst, clusterMST) 
  }
  
  if (8 %in% pipeline){
    mst_visualize(colormst, bordermst, clusterMST, ig_mst$ig, ig_mst$cbg, ig_mst$cbor, ig_mst$coords)
    mst_with_clustering_colors(colormst, bordermst, clusterColor, centralColor, clusterSelect2, clusterMST)
  }
  
  #Centrality 
  if (9 %in% pipeline){
    calculate_centrality(fastCButton, centralSelect, clusterCentral, centralColor)
  }
  
  #Label propagation - classification
  if (10 %in% pipeline){
    w <- as.matrix(values$sim)
    #all_ids <- values$indexes
    known.label <- which(!is.na(data$V2))
    y <- data$V1[known.label]
    
    if (9 %in% pipeline){
      c_s <- centralValues$summary$Mean/max(centralValues$summary$Mean)
      c_s <- 0 # not use this
    }else{
      c_s <- 0
    }
    
    f1<-sslLabelProp(w,y,known.label,1000,c_s)
    t <- f1$t
    f1 <- f1$f
    
    #False 
    a1 <- length(which(abs(data$V1[-known.label]-f1[-known.label])<1))/length(as.vector(data$V1[-known.label])) *100
    print(a1)
    print(paste0("Num of elements with diff from actual class < 1 : ", length(which(abs(data$V1[-known.label]-f1[-known.label])<1)) ))
    print(paste0("Num of truly classified elements: ", length(which(round(f1[-known.label])==as.vector(data$V1[-known.label])))))

    a <- NA
    tryCatch({
      a <- accuracy(x=as.vector(round(f1[-known.label])),y=as.vector(data$V1[-known.label]))$PCC
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
#print("x-y=")
#print( (as.vector((f1[-known.label]))-as.vector(data$V1[-known.label]))) 
print(paste0("Accuracy: ", a))
    classification_eval_matrix <- as.data.frame(matrix(0, nrow = 1, ncol = 5))
    colnames(classification_eval_matrix) <- c("Algorithm", "N_actual", "Accuracy", "Accuracy_with_margin_1", "Time")
    classification_eval_matrix$Algorithm[1] <- "label.propagation"
    classification_eval_matrix$N_actual[1] <- length(unique(y))
    classification_eval_matrix$Accuracy[1] <- a
    classification_eval_matrix$Accuracy_with_margin_1[1] <- a1
    classification_eval_matrix$Time[1] <- t

    if(!file.exists(paste0(output_folder,"/classification"))){ 
      dir.create(paste0(output_folder,"/classification"))
    }
    write.table(classification_eval_matrix, file = paste0(output_folder, "/classification", "/", "classification_", 
                                                          with_edges, "_", graph.type,"_", thr, "_", iter, ".txt"), row.names = F)
    
  }
 
  #f1<-sslLabelProp(w,y,known.label,graph.type="enn",epsilon = 0.5)
  #f2<-sslLabelProp(w,y,known.label,graph.type="knn",k =10)
  #f3<-sslLabelProp(w,y,known.label,graph.type="tanh",alpha1=-2,alpha2=1)
  #f4<-sslLabelProp(w,y,known.label,graph.type="exp",alpha = 1)
  
  ##Performs MST-kNN clustering using euclidean distance
  #results <- mst.knn(w)
  ## Visualizes the clustering solution
  #library("igraph")
  
  #png(paste0(output_folder,"/","knn.mst.png"))
  #plot(results$network, vertex.size=8,
  #     vertex.color=igraph::clusters(results$network)$membership,
  #     layout=igraph::layout.fruchterman.reingold(results$network, niter=10000),
  #     main=paste("MST-kNN \n Clustering solution \n k=",results$k,sep="" ))
  #dev.off()
  
}
