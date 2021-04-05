#!/usr/bin/env Rscript
# Script for Pooling cells having ERB variants
# @Arpan Kumar Basak

# Script
rm(list = ls())
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
path <- as.character(getwd())
setwd(path)

# For colour donot use thresholding, instead treat the red channel as potential cellwall eroding the lines
message(paste("---", "== Pooling features =="," Job Started at ", as.character(Sys.time()), "---", sep = ""))

# Load Requisites
pkgs <- c("EBImage", "tidyverse", "parallel")
lapply(pkgs, require, character.only = T)
options(warn = 1, mc.cores = 8)

# Measures not to fry R
options(set.cores = 8)
param <- "./script/parameters.R"
source(param)

cell_path <- paste("./images/segmented_images/stacked/cells/")

load(paste0(data, "Clustered_feature_image_data.RData"))

# Fetch micro-features
fls_cell <- list.files(path = cell_path,
                  full.names = F, 
                  include.dirs = F, 
                  ignore.case = T, 
                  recursive = F, 
                  pattern = ".png")

# Filter out garbage
id_cell <- (fls_cell %in% paste0("cell_", imat.clustered$FeatureID, ".png"))
fls_cell <- fls_cell[id_cell]

# Represent the images for each clusters and each genotype
mclapply(unique(as.character(imat.clustered$k_cluster)), function(x){

    # x <- "cluster60"
    clust_path <- paste(path, "/images/clustered_images/", x, sep = "")
    if(!dir.exists(clust_path)){
      dir.create(path = clust_path, recursive = TRUE)
    }
    temp <- imat.clustered %>% filter(k_cluster == x)
    idx <- as.character(temp$cell_path)
    idx <- fls_cell[fls_cell %in% idx]
    temp$file_path <- paste0(cell_path, idx[match(temp$cell_path, idx)])

    if(!is.na(length(idx))){

      # Make stacked images
      message(paste0("Files in ", x, " : "))
      message(length(unique(temp$file_path)))
      for(j in 1:length(unique(temp$file_path))){

          i <- unique(temp$file_path)[j]
          img <- readImage(files = i, type = "png")
          tx <- str_split(i, pattern = "/", simplify = T)[,6]

           # Write stacked images
           writeImage(img, paste(clust_path, "/", x, "_", tx,
                                            sep = ""), 
           quality = 100)

      }
      message(paste0("Done saving ", x, " : "))
      
      }else {
        
        message("NONE")
      }
    
   
}, mc.cores = 4)

# Prepare the datatable in html format
imat.clustered$cell_path_clustered <- img_uri(imat.clustered$cell_path)
tab <- DT::datatable(imat.clustered, 
  filter = 'top', 
  options = list(pageLength = 5, autoWidth = TRUE))

htmlwidgets::saveWidget(tab, "./output/clustered_cells.html")


message("DONE!!")
sessionInfo()