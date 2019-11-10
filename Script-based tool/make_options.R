#!/usr/bin/env Rscript

# Wrapper for the run_graph_tool_without_ui script. It lives as an R shell script in 
# order to be seamlessly integrated to the general variant calling pipeline 
# which uses only system interactions.
# To run the script using the command line using the default input arguments: Rscript --vanilla make_options.R
# e.g. Rscript --vanilla make_options.R --pipeline 4
# Help: Rscript --vanilla make_options.R --help 

# Author: Maria Kotouza

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(
    opt_str=c("-w","--dataset_name"),
    action="store",
    default="tr23.wc",
    help=paste0(
      "The directory where the input file is located."
    )
  ),
  make_option(
    opt_str=c("-a","--datapath"),
    action="store",
    default="data_discr/tr23.wc/udata_classes_grouped.txt",
    help=paste0(
      "The directory where the input file is located."
    )
  ),
  make_option(
    opt_str=c("-b","--sep"),
    action="store",
    default="\t",
    help=paste0(
      "The files of the IMGT output that will be used through the analysis.\n",
      "Use comma to separate the list of files (default %default)"
    )
  ),
  make_option(
    opt_str=c("-c","--group_by_clusterId"),
    action="store",
    default="F",
    help="Group by sequence-cluster id: 'T' for True or 'F' for False."
  ),
  make_option(
    opt_str=c("-v","--delete_edges_of_diff_classes"),
    action="store",
    default="T",
    help="Delete edges that connect nodes of different classes: 'T' for True or 'F' for False."
  ),
  make_option(
    opt_str=c("-d","--seqSelect"),
    action="store",
    default="AA.JUNCTION",
    help="Select the column of the sequence to calculate the distance matrix."
  ),
  make_option(
    opt_str=c("-e","--simSelect"),
    action="store",
    default="lv",
    help=paste0(
      "Select the similarity metric:\n",
      "osa for OSA \n",
      "lv for Levenshtein \n",
      "dl for Full DL \n", 
      "hamming for Hamming \n",
      "lcs for LCS \n",
      "qgram for Q-Gram \n",
      "cosine for Cosine \n",
      "jaccard for Jaccard \n",
      "jw for Jaro Winklar \n"
    )
  ),
  make_option(
    opt_str=c("-s","--graph.type"),
    action="store",
    default="knn:10",
    help=paste0(
      "Select the sparsification methods to be used. Select between 'knn', 'enn', 'tanh', 'exp'. \n",
      "For knn specify the number of neighbors using ':' as follows knn:10 \n",
      "For mknn specify the number of neighbors using ':' as follows mknn:10 \n",
      "For enn specify the e-neighbor threshold using ':' as follows enn:0.5. The threshold takes values from [0,1] - 1 holds all the connections, whereas 0 deletes all the connections \n",
      "default %default"
    )
  ),
  make_option(
    opt_str=c("-n","--pipeline"),
    action="store",
    default="4",
    help=paste0(
      "Pipeline options:\n",
      "1. Filtering \n",
      "2. Visualize Graph \n",
      "3. Select component \n",
      "4. Clustering \n",
      "5. Visualize clustered graph \n",
      "6. Compare clustering results between two methods \n",
      "7. MST construction \n",
      "8. Visualize MST \n",
      "9. Centrality Computation \n",
      "10. Classification \n",
      "Use comma to separate the list of options (default %default)"
    )
  ),
  make_option(
    opt_str=c("-g","--componentSelect"),
    action="store",
    default="1",
    help=paste0(
      "Select the component id to be used for the rest of the analysis. \n"
    )
  ),
  make_option(
    opt_str=c("-i","--excludeButton"),
    action="store",
    default="F",
    help=paste0(
      "Filter out rows from the dataset. Select Yes or No.\n"
    )
  ),
  make_option(
    opt_str=c("-j","--includeButton"),
    action="store",
    default="F",
    help=paste0(
      "Filter in rows from the dataset. Select Yes or No.\n"
    )
  ),
  make_option(
    opt_str=c("-k","--select1"),
    action="store",
    default="AA.JUNCTION",
    help=paste0(
      "The attribute that you want to apply a filter on.\n"
    )
  ),
  make_option(
    opt_str=c("-l","--slider1"),
    action="store",
    default="1",
    help=paste0(
      "This is used in case that the attribute select1 of the filtering process is numeric. Insert the value of the attribute.\n"
    )
  ),
  make_option(
    opt_str=c("-m","--text1"),
    action="store",
    default="",
    help=paste0(
      "This is used in case that the attribute select1 of the filtering process is categorical. Insert the value of the attribute.\n"
    )
  ),
  make_option(
    opt_str=c("-o","--mstAlgoSelect"),
    action="store",
    default="Prim",
    help=paste0(
      "The MST algorithm\n",
      "Select 'Prim' or 'Kruskal'\n",
      "default %default"
    )
  ),
  make_option(
    opt_str=c("-p","--clusterMST_param"),
    action="store",
    default="T;V1;V1",
    help=paste0(
      "Select 'T' for True - color according to clustering results, or 'F'.\n",
      "Select the attribute to use for node coloring. Use Default for not coloring.\n",
      "Select the attribute to use for boarder-node coloring. Use Default for not coloring.\n",
      "Seperate the options using semicolon ;",
      "default %default"
    ) 
  ),
  make_option(
    opt_str=c("-q","--centralities"),
    action="store",
    default="F;Shortest.Paths.Betweenness.Centrality;F;AA.JUNCTION",
    help=paste0(
      "Select 'T' for True - Only major centralities metrics, or 'F' - compute all centrality metrics.\n",
      "Select a centrality type between 'Degree.Centrality', 'Shortest.Paths.Betweenness.Centrality', 'Average.Distance', 'eigenvector.centralities'.\n",
      "Select 'T' for True - color according to clustering results, or 'F'. \n",
      "Select the attribute to use for node coloring. Use Default for not coloring.\n",
      "Seperate the options using semicolon ;",
      "default %default"
    )
  ),
  make_option(
    opt_str=c("-r","--clustering"),
    action="store",
    default="louvain;fast_greedy;label_propagation;leading_eigenvalue;walktrap",
    help=paste0(
      "Select the different clustering methods to be used. Select between 'louvain', 'fast_greedy', 'label_propagation', 'leading_eigenvalue', 'walktrap', 'edge_betweenness', 'hierarchical'. \n",
      "Select 'reads' of 'theshold' for clonotypes, the number of reads or the threshold percentage accordingly, and whether you want to take gene into account\n",
      "Seperate the methods using semicolon ;\n",
      "default %default"
    )
  ),
  make_option(
    opt_str=c("-x","--clusterColor"),
    action="store",
    default="Default",
    help=paste0(
      "Select the attribute to use for node coloring. Use Default for coloring according to the clustering results.\n",
      "default %default"
    )
  ),
  make_option(
    opt_str=c("-t","--plotHeat"),
    action="store",
    default="T",
    help=paste0(
      "\n",
      "default %default"
    )
  ),
  make_option(
    opt_str=c("-u","--heatSelect"),
    action="store",
    default="0",
    help=paste0(
      "\n",
      "default %default"
    )
  ),
  make_option(
    opt_str=c("-y","--percentage_labeled_data"),
    action="store",
    default="0.30",
    help=paste0(
      "Select the percentage of labeled data. Choose between [0,1] \n",
      "default %default"
    )
  ),
  make_option(
    opt_str=c("-i","--iter"),
    action="store",
    default="1",
    help=paste0(
      "Iteration for cross-validation \n",
      "default %default"
    )
  )
);

#other options that have not been used: -f


opt <- parse_args(OptionParser(option_list=option_list));


source("run_graph_tool_without_ui.R")
run_graph_tool(dataset_name = opt$dataset_name, datapath=opt$datapath, sep=sep, header=T, pipeline=opt$pipeline, 
               group_by_clusterId=as.logical(opt$group_by_clusterId), seqSelect=opt$seqSelect, graph.type=strsplit(opt$graph.type, ";")[[1]],
               simSelect=opt$simSelect, componentSelect=as.numeric(opt$componentSelect), 
               excludeButton=as.logical(opt$excludeButton), includeButton=as.logical(opt$includeButton), select1=opt$select1,
               slider1=as.numeric(opt$slider1), text1=opt$text1, mstAlgoSelect=opt$mstAlgoSelect, 
               delete_edges_of_diff_classes=as.logical(opt$delete_edges_of_diff_classes),
               clusterMST=as.logical(strsplit(opt$clusterMST_param, ";")[[1]][1]), colormst=strsplit(opt$clusterMST_param, ";")[[1]][2],
               bordermst=strsplit(opt$clusterMST_param, ";")[[1]][3],
               fastCButton=as.logical(strsplit(opt$centralities, ";")[[1]][1]), centralSelect=strsplit(opt$centralities, ";")[[1]][2], 
               clusterCentral=as.logical(strsplit(opt$centralities, ";")[[1]][3]), centralColor=strsplit(opt$centralities, ";")[[1]][4],
               clustering_algo_Select=strsplit(opt$clustering, ";")[[1]],
               clusterColor=opt$clusterColor, plotHeat=as.logical(opt$plotHeat), heatSelect=as.numeric(opt$heatSelect),
               percentage_labeled_data=as.numeric(opt$percentage_labeled_data), iter=as.numeric(opt$iter))
# Print help
# print_help( OptionParser(option_list=option_list))

