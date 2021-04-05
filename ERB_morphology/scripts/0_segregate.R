#!/usr/bin/env Rscript

rm(list = ls())
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
path <- as.character(getwd())

setwd(path)

# Run this script for annotating the grayscale images only
# For colour donot use thresholding, instead treat the red channel as potential cellwall eroding the lines
message(paste("---", "==Segregation=="," Job Started at ", as.character(Sys.time()), "---", sep = ""))
# Load Requisites
pkgs <- c("EBImage", "MaxContrastProjection", "tidyverse", "parallel", "stringr")
lapply(pkgs, require, character.only = T)

# Measures not to fry R
options(set.cores = detectCores()/2)
source("./script/parameters.R")

# Path arguments to be taken from the user command
image_stock = as.character(args[1]) # Path to the image directory
param = as.character(args[2]) # Path to the parameter

# List the directories
dir_list <- c(image, 
  dyn_image,
  cell_dyn_image,
  live_stock, 
  cell_image, 
  feature_image, 
  live_stock, 
  output, 
  feature_mat,
  stats, 
  figs
  )

# Directory check and create directories
mclapply(dir_list, function(x){

  if(!dir.exists(x)){

  message("Directory created !!")

  dir.create(path = x, recursive = TRUE)

} else{

  message("Directory exists !!")
}
}, mc.cores = 4)


# List files
df_param <- read.table(
  file = "./data/global_parameters.txt", 
  stringsAsFactors = FALSE, 
  header = TRUE, 
  sep = "\t"
) %>% 
filter(mode %in% c("segregate", "cores")) %>%
data.frame(., stringsAsFactors = FALSE, row.names = .$parameter) %>% 
select(default) %>%
t() %>%
data.frame(., stringsAsFactors = FALSE)

# List files
fls <- list.files(path = image_stock,
                  full.names = T, 
                  include.dirs = T, 
                  ignore.case = T, 
                  recursive = T, 
                  pattern = ".tif")

head(fls)
message("...")
tail(fls)

# Subset parameter 2 and compute features
# idx <- str_detect(fls, "_param1_")
length(fls)
fls <- as.list(fls)

id <- str_split(fls, pattern = "_z", simplify = T)[,1]
id <- str_replace(id, pattern = image_stock, replacement = "")
id <- unique(str_split(id, 
                       pattern = "/", simplify = T)[,2])

ch_red <- as.numeric(df_param$channel_red)
ch_green <- as.numeric(df_param$channel_green)
ch_blue <- as.numeric(df_param$channel_blue)
dim_img <- as.numeric(df_param$dimension)
mcores <- as.numeric(df_param$mc.cores)

mclapply(1:length(id), function(i) {
    
    message(paste("Processing file ", id[i], ""))

    set.seed(1)
    ix <- grep(fls, pattern = id[i], value = TRUE)
    
    # Load Image - set for ERBs and Cells - Fetch a concat image for feature and cell
    img <- readImage(unlist(ix), type = "tiff")
    
    if(dim(img)[2] > dim_img) {
        message(paste0("Image Dimension: ", dim(img)[1], " x ", dim(img)[2]))
        img <-img[,1:dim_img,,]
    }
     else {
        message(paste0("Image Dimension: ", dim(img)[1], " x ", dim(img)[2]))
        feature <- img[,,ch_red,]
        cell <- img[,,ch_green,]
        colorMode(img) <- Grayscale
    }
    
    # Conduct MaxContrastProjection
    cell <- contrastProjection(imageStack = cell, 
                                 w_x = as.numeric(df_param$cell_window_x),
                                 w_y = as.numeric(df_param$cell_window_y), 
                                 brushShape = as.character(df_param$cell_brush_shape), 
                                 smoothing = as.numeric(df_param$cell_smoothing))
        
    feature <- contrastProjection(imageStack = feature, 
                                 w_x = as.numeric(df_param$feature_window_x),
                                 w_y = as.numeric(df_param$feature_window_y), 
                                 brushShape = as.character(df_param$feature_brush_shape),
                                 smoothing = as.numeric(df_param$feature_smoothing))
    
    # Output Image
    img <- rgbImage(red = cell, green = feature)[,,-3]
    colorMode(img) <- Grayscale
    img  <- normalize(img)
         
    message(paste("\nKnitting image ", id[i], sep = ""))
    writeImage(x = img, 
               files =  paste(live_stock, id[i], ".tif", sep = ""), 
               quality = 100, 
               type = "tiff", 
               bits.per.sample = 8)
    
    message("DONE!!")
    
}, mc.cores = mcores)

message(paste("---", "==Segregation=="," Job ended at ", as.character(Sys.time()), "---", sep = ""))

# END OF SCRIPT!
sessionInfo()
