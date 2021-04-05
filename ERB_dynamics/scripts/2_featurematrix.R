#!/usr/bin/env Rscript

# Script for Cell image feature segmentation
# @Arpan Kumar Basak
rm(list = ls())
options(warn = 1)
param <- "./script/parameters.R"
path <- as.character(getwd())

# Run this script for annotating the grayscale images only
# For colour donot use thresholding, instead treat the red channel as potential cellwall eroding the lines
message(paste("---", "==Building Feature Matrix=="," Job Started at ", as.character(Sys.time()), "---", sep = ""))

# Load Requisites
pkgs <- c("EBImage", "tidyverse", "parallel")
lapply(pkgs, require, character.only = T)

# Measures not to fry R
options(set.cores = 8)
source(param)
setwd(path)
nc <- 8 # Partition cores

# List files
fls <- list.files(path = paste(output, "feature_matrix/", 
                               sep = ""), 
                  full.names = T, include.dirs = F)
metadata <- data.frame()
imat <- data.frame()
dat <- data.frame()

# Loop to list out the feature data from images processed
res <- mclapply(fls, function(x){
    
    message(paste("Reading feature data of ", x, sep = ""))
    temp <- read.delim2(x, header = T, sep = "\t", as.is = T)
    if(nrow(temp) > 1){
        message("Done!!")
    } else {
        message("None features detected in the segmentation process!!")
    }
    return(temp)
}, mc.cores = 24)

# Use names for each matrix from filename
id <- str_replace(string = fls, 
                  pattern = paste(output, "feature_matrix/", sep = ""), 
                  replacement = "")
id <- str_replace_all(string = id, pattern = "/", replacement = "")
id <- str_replace_all(string = id, pattern = ".txt", replacement = "")
names(res) <- id


# Make two matrices for cell features and object features
idx <- grep("^cell", names(res), ignore.case = T)
imat_c <- do.call(rbind, res[idx])

idx <- grep("^feature", names(res), ignore.case = T)
imat <- do.call(rbind, res[idx])

rm(list = "res") # Remove list object to save RAM


message("Glimpse...")
head(imat)
tail(imat)
message("Dimension...")
dim(imat)

message("Saving output...")

save(list = c("imat"), file = paste(data, "Obj_image_features.RData", sep = ""))

# END OF SCRIPT
sessionInfo()
