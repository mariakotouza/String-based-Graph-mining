#!/usr/bin/env Rscript

# Wrapper for the run_discretization script. It lives as an R shell script in 
# order to be seamlessly integrated to the general variant calling pipeline 
# which uses only system interactions.
# To run the script using the command line using the default input arguments: Rscript --vanilla make_options_discretization.R

# Author: Maria Kotouza

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(
    opt_str=c("-a","--datapath"),
    action="store",
    default="Data/tr23.wc/20topics",
    help=paste0(
      "The directory where the dataset is located."
    )
  ),
  make_option(
    opt_str=c("-b","--dataset_name"),
    action="store",
    default="tr23.wc",
    help="The name of the dataset."
  ),
  make_option(
    opt_str=c("-c","--numeric_file_name"),
    action="store",
    default="model-final.theta",
    help=paste0(
      "The name of the file with the numeric values seperated by space. Each raw contains the observations of an item.\n",
      "(default %default)"
    )
  ),
  make_option(
    opt_str=c("-d","--ids_file_name"),
    action="store",
    default="docIDs.dat",
    help=paste0(
      "The ids of the items. Each raw contains an item id.\n",
      "(default %default)"
    )
  ),

  make_option(
    opt_str=c("-e","--num_of_bins"),
    action="store",
    default="10",
    help=paste0(
      "The number of bins selected for discretization (default %default) \n"
    )
  ),
  make_option(
    opt_str=c("-f","--num_of_bins_per_group"),
    action="store",
    default="2",
    help=paste0(
      "The number of consecutive bins that will be assigned to a group for the Similarity-based algorithm (default %default)"
    )
  ),
  make_option(
    opt_str=c("-g","--num_of_topics"),
    action="store",
    default="20",
    help="The number of topics."
  )
);

opt <- parse_args(OptionParser(option_list=option_list));

packages <- c("plyr","dplyr","data.table","stringr","tidyr","entropy","ggplot2","ggseqlogo","gridExtra","cluster","stringdist","igraph","data.tree", "networkD3","xtable","tictoc","parallel")
lapply(packages, library, character.only = TRUE)

source("functions.R")
source("run_discretization.R")
run_discretization(datapath=opt$datapath, dataset_name=opt$dataset_name, numeric_file_name=opt$numeric_file_name, 
                   ids_file_name=opt$ids_file_name,
                   num_of_bins=as.numeric(opt$num_of_bins), num_of_bins_per_group=as.numeric(opt$num_of_bins_per_group),
                   num_of_topics=as.numeric(opt$num_of_topics))

# Print help
# print_help( OptionParser(option_list=option_list))

