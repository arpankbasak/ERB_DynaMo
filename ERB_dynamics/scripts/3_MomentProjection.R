#!/usr/bin/env Rscript

rm(list = ls())
options(warn = 1)
# args <- commandArgs(trailingOnly = TRUE)
path <- as.character(getwd())
setwd(path)

# Run this script for annotating the grayscale images only
# For colour donot use thresholding, instead treat the red channel as potential cellwall eroding the lines
message(paste("---", "==MomentProjection=="," Job Started at ", as.character(Sys.time()), "---", sep = ""))
# Load Requisites
pkgs <- c("EBImage", "tidyverse", "parallel", "stringr", "imager")
lapply(pkgs, require, character.only = T)

# Measures not to fry R
options(set.cores = detectCores()/2)


# Set path
image_frames = "./images/merged_images"
cell_images = "./images/segmented_images/stacked/cells/"
dyn = "./images/dynamics/"
dyn_cell = "./images/dynamics/cells/"

# List files
fls <- list.files(path = image_frames,
                  full.names = T, 
                  include.dirs = T, 
                  ignore.case = T, 
                  recursive = T, 
                  pattern = ".png")

head(fls)
message("...")
tail(fls)

# Subset parameter 2 and compute features
# fls <- as.list(fls[!idx])

id <- str_split(fls, pattern = "_z", simplify = T)[,1]
id <- str_replace(id, pattern = image_frames, replacement = "")
# id <- str_replace(id, pattern = eval_stock, replacement = "")
id <- unique(str_replace_all(id, 
                       pattern = "/", ""))

# Make video of the time series snaps
system(paste0("rm -rf ", dyn, "*"))
mclapply(1:length(id), function(i) {
    
    temp <- id[i]
    ix <- grep(fls, pattern = temp, value = TRUE, ignore.case = TRUE)
    av::av_encode_video(ix, paste0(dyn, temp, ".mp4"), 
    	framerate = 2)
   
    
}, mc.cores = 8)

message(paste("---", "==Dynamics are projected=="," Job ended at ", as.character(Sys.time()), "---", sep = ""))

# END OF SCRIPT!
sessionInfo()
