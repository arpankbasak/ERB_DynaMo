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

time_points <- as.numeric(df_param$time_points)
tim <- paste0("_t", 0:(time_points-1))

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
# fls <- as.list(fls[!idx])

id <- str_split(fls, pattern = "_z", simplify = T)[,1]
id <- str_replace(id, pattern = image_stock, replacement = "")
# id <- str_replace(id, pattern = eval_stock, replacement = "")
id <- unique(str_split(id, 
                       pattern = "/", simplify = T)[,3])

# idx <- sample(id, 5)

ch_red <- as.numeric(df_param$channel_red)
ch_green <- as.numeric(df_param$channel_green)
ch_blue <- as.numeric(df_param$channel_blue)
dim_img <- as.numeric(df_param$dimension)
mcores <- as.numeric(df_param$mc.cores)

mclapply(1:length(id), function(i) {
    
    # TODO::Add a conditional statement to check the channel wether gray or rgb
    # TODO::Include Messages at each step of analysis for log and progress
    message(paste("Processing file ", id[i], ""))
    # x <- (fls[1])[11]
    # i <- 11
    # j <- 4
    lapply(1:length(tim), function(j){

        set.seed(1)
        flt <- grep(fls, pattern = tim[j], value = TRUE)
        ix <- grep(flt, pattern = id[i], value = TRUE)
        
        # Load Image - set for ERBs and Cells - Fetch a concat image for feature and cell
        img <- readImage(unlist(ix), type = "tiff")
        dp <- dim(img)[4]
        
        if(dim(img)[2] > dim_img) {
            message(paste0("Image Dimension: ", dim(img)[1], " x ", dim(img)[2]))
            img <-img[,1:dim_img,,]
        }else {
            message(paste0("Image Dimension: ", dim(img)[1], " x ", dim(img)[2]))
            feature <- img[,,ch_green,]
            cell <- img[,,ch_red,]
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
             
        message(paste("\nKnitting image ", id[i], tim[j], sep = ""))
        writeImage(x = img, 
                   files =  paste(live_stock, id[i], "_z", dp, tim[j], ".tif", sep = ""), 
                   quality = 100, 
                   type = "tiff", 
                   bits.per.sample = 8)
    
    message("DONE!!")
    
        })
    
}, mc.cores = mcores)

# Find the N time and integrate in B channel for moment dynamics
fls_live <- list.files(path = live_stock,
                  full.names = T, 
                  include.dirs = T, 
                  ignore.case = T, 
                  recursive = T, 
                  pattern = ".tif")

id <- str_split(fls_live, pattern = "_t", simplify = T)[,1]
id <- str_split(fls_live, pattern = "_z", simplify = T)[,1]
id <- str_replace(id, pattern = live_stock, replacement = "")
# id <- unique(str_split(id, 
#                        pattern = "/", simplify = T)[,3])

mclapply(1:length(unique(id)), function(i){

  temp <- unique(id)[i]

  idx <- str_detect(fls_live, temp)
  all_cards <- fls_live[idx]

  id0 <- str_detect(all_cards, "_t0")
  ace <- all_cards[id0]

  # fetch the ace 
  img <- readImage(ace, type = "tiff")
  initial_position <- img[,,ch_green]

  lapply(1:length(all_cards), function(j){

    slice <- all_cards[j]
    img <- readImage(slice, type = "tiff")
    img <- rgbImage(red = img[,,ch_red], 
      green = img[,,ch_green], 
      blue = initial_position)[,,]
    colorMode(img) <- Grayscale
    writeImage(x = img, 
                   files =  paste(slice, sep = ""), 
                   quality = 100, 
                   type = "tiff", 
                   bits.per.sample = 8)

  })


}, mc.cores = mcores)


message(paste("---", "==Segregation=="," Job ended at ", as.character(Sys.time()), "---", sep = ""))

# END OF SCRIPT!
sessionInfo()
