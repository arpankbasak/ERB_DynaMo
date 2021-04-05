#!/usr/bin/env Rscript

# Script
rm(list = ls())
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
path <- as.character(getwd())
setwd(path)

# Initiate a jupyter instance for evaluating the dynamics

# Run this script for annotating the grayscale images only
message(paste("---", "==Segmentation=="," Job Started at ", as.character(Sys.time()), "---", sep = ""))

# Load Requisites
pkgs <- c("EBImage", "tidyverse", "parallel")
lapply(pkgs, require, character.only = T)

# Measures not to fry R
options(set.cores = 8)
param <- "./script/parameters.R"
source(param)

# param = as.character(args[1]) # Path to the parameter
param = "./data/global_parameters.txt"
df_param <- read.table(
  file = param, 
  stringsAsFactors = FALSE, 
  header = TRUE, 
  sep = "\t"
) %>% 
filter(mode %in% c("segmentation", "cores")) %>%
data.frame(., stringsAsFactors = FALSE, row.names = .$parameter) %>% 
select(default) %>%
t() %>%
data.frame(., stringsAsFactors = FALSE)

mcores <- as.numeric(df_param$mc.cores)


# List files
fls <- list.files(path = live_stock, 
                  full.names = T, 
                  include.dirs = T, 
                  ignore.case = T, 
                  recursive = T, 
                  pattern = ".tif")

message(paste("No. of images: ", length(fls), sep = ""))

res <- list()
res <- mclapply(1:length(fls), function(ind) {
    
    dat <- list()
    out = output
    x <- fls[ind]
    set.seed(1)
    
    # Load Image - set for ERBs
    message(paste("Processing file ", x, ""))
    img <- readImage(x, type = "tiff")
    if(validObject(img)){


        
        brush <- makeBrush(size=as.numeric(df_param$kernal_brush_size), 
            shape = as.character(df_param$kernal_brush_shape), 
            sigma=as.numeric(df_param$kernal_brush_sigma))
        feature <- filter2(img[,,2], filter = brush)
        # Morphology operation for features and cell

        message("Processing cellular features")
        feature <- feature * (feature > otsu(feature))
        # feature <- feature * (feature > quantile(feature, as.numeric(df_param$feature_q_threshold)))
        feature <- thresh(feature, 
            w = as.numeric(df_param$feature_opening_w), 
            h = as.numeric(df_param$feature_opening_h), 
            offset = as.numeric(df_param$feature_opening_offset))
        feature <- opening(feature, 
            kern = makeBrush(as.numeric(df_param$feature_opening_brush_width), 
            shape= as.character(df_param$feature_opening_brush_shape)))
        eseed <- bwlabel(feature)
        #feature <- fillHull(thresh(temp, w = 20, h = 20, offset = 0.00005))
        emask <- propagate(img[,,2], eseed, mask=feature)

        # Extract cell features
        message(paste("Processing Cell Wall like features", "..."))

        # Adaptive thresholding for vacuolar space
        vac.f <- filter2(img[,,1], filter = brush)
        vac <- vac.f
        vac[!(vac.f > 0 & vac.f < quantile(vac.f, 0.5))] <- -1
        vac <- thresh(vac, w = 20, h = 20, offset = 0.000001)
        vac <- opening(vac, kern=makeBrush(3,shape="disc"))
        vac.o <- opening(!(vac.f > 0 & vac.f < quantile(vac.f, 0.65)))
        vac[vac.o>vac] <- vac[vac.o>vac]
        vac <- fillHull(vac)
        vseed <- bwlabel(vac)
        vmask <- propagate(img[,,1], vseed, mask=vac)

        # Adaptive thresholding for cell and ER-net space
        vac.f <- filter2(img[,,1], filter = brush)
        vac <- vac.f
        vac[!(vac.f > 0 & vac.f < quantile(vac.f, as.numeric(df_param$void_brush_q50)))] <- -1
        vac <- thresh(vac, 
            w = as.numeric(df_param$void_opening_w), 
            h = as.numeric(df_param$void_opening_h), 
            offset = as.numeric(df_param$void_opening_offset))
        vac <- opening(vac, 
            kern=makeBrush(as.numeric(df_param$void_opening_brush_width),
            shape= as.character(df_param$void_opening_brush_shape)))
        vac.o <- opening(!(vac.f > 0 & vac.f < quantile(vac.f, as.numeric(df_param$void_brush_qrange))))
        vac[vac.o>vac] <- vac[vac.o>vac]
        vac <- fillHull(vac)
        vseed <- bwlabel(vac)
        vmask <- propagate(img[,,1], vseed, mask=vac)

        # Adaptive thresholding for cell and ER-net space
        cell <- filter2(img[,,1], filter = brush)
        cell <- cell * (cell < quantile(cell, 0.95))
        cell <- fillHull(thresh(cell, 
            w = as.numeric(df_param$cell_hull_h), 
            h = as.numeric(df_param$cell_hull_w), 
            offset = as.numeric(df_param$cell_opening_offset)))
        cell <- opening(cell, kern=makeBrush(as.numeric(df_param$cell_opening_brush_width),
            shape=as.character(df_param$cell_opening_brush_shape)))
        cell.o <- opening(cell < quantile(cell, as.numeric(df_param$cell_opening_qlimit)))
        cell[cell.o>cell] <- cell.o[cell.o>cell]
        cell[vac > cell] <- vac[vac > cell]
        cell[feature > cell] <- feature[feature > cell]
        cmask <- propagate(img[,,1], vmask, lambda=as.numeric(df_param$cell_propagate_lambda), mask=cell)

        # Segmentation
        seg_obj <- paintObjects(vac, img[,,1], thick = T, col = "#ffff31")
        seg_obj <- paintObjects(cell, seg_obj, thick = F, col = "#ffff10")
        seg_obj <- paintObjects(feature, seg_obj, thick = F, col = "#ffff00")

        merge_obj <- rgbImage(red = img[,,1] * cell, 
                             green = img[,,2] * feature, 
                             blue = img[,,1] * vac)
        
        # Merge and paint the cell-like features and stack them
        cell_stack_list <- EBImage::stackObjects(cmask, merge_obj)
        # display(cell_stack_list[,,-id_fpc], all = TRUE)

        # Get cell-features
        a <- as.data.frame(computeFeatures.shape(cmask, merge_obj))
        b <- as.data.frame(computeFeatures.moment(cmask, merge_obj))
        c <- as.data.frame(computeFeatures.haralick(cmask, merge_obj))
        d <- as.data.frame(computeFeatures.basic(cmask, merge_obj))

        cell_stack_list <- EBImage::stackObjects(cmask, merge_obj)
        # display(cell_stack_list[,,-id_fpc], all = TRUE)

        # Get cell-features
        a <- as.data.frame(computeFeatures.shape(cmask, merge_obj))
        b <- as.data.frame(computeFeatures.moment(cmask, merge_obj))
        c <- as.data.frame(computeFeatures.haralick(cmask, merge_obj))
        d <- as.data.frame(computeFeatures.basic(cmask, merge_obj))

        cdat <- cbind(a, b, c, d)
        flag = FALSE
        if(nrow(cdat) != 0){
            flag = TRUE
        }

        if(flag == TRUE){
            print(dim(cdat))
            print(head(cdat))
            message("...")
            print(nrow(cdat))
        }else {
            message("None cells detected; change parameter if required!!")
            }


        message("==== Take a glimpse ====")

        # Store SampleID
        sample_id <- str_replace(x, pattern = ".tif", replacement = "")
        sample_id <- str_replace(sample_id, pattern = as.character(live_stock), replacement = "")
        sample_id <- str_replace(sample_id, pattern = "/", replacement = "")
        writeImage(merge_obj, paste("./images/merged_images/", 
                                         sample_id, ".png", 
                                          sep = ""), quality = 100)
        # Store meta information
#         etemp = NULL
        if(nrow(cdat) != 0){

            # Total cell like features
            cdat$total_counts <- nrow(cdat)

            # Obtain the feature of confidence
            id_fpc <- which(cdat[,1] > as.numeric(df_param$cell_rarefry_surface_area)) # Filter garbage or rarefy
            cdat$cell_feature <- 0
            cdat$cell_feature[id_fpc] <- 1
            if(length(id_fpc) != 0){
                
                cdat <- cdat %>% filter(cell_feature == 1) %>% add_column(SampleID = sample_id, .before = 1)
                cdat <- cdat %>% add_column(FeatureID = paste(sample_id, row.names(cdat), sep = "_"),
                                            .before = 1)
                cell_stack_list <- cell_stack_list[,,,id_fpc]
                message("Filtered potential cells:")
                glimpse(cdat)
                # Potential Cells
                wx <- colorLabels(cmask)
                writeImage(wx, paste("./images/segmented_images/cell_binary_segmentation-", 
                                     sample_id, ".png", sep = ""), quality = 100)

                # Merged and segmented
                wx <- merge_obj
                writeImage(wx, paste("./images/segmented_images/merged_segmentation-", 
                                     sample_id, ".png", 
                                     sep = ""), quality = 100)

                # Count number of ERB signals detected
                cdat$nfeatures <- 0
                edat <- list()
                for(cell_idx in 1:nrow(cdat)){

                        # Segment ERbodies
        #                 cell_idx <- 13
                        mask <- cell_stack_list[,,2, cell_idx]
                        colorMode(mask) <- Grayscale
                        mask <- mask * (mask > otsu(mask))
#                         mask <- mask * (mask > quantile(mask, 0.9))
                        mask <- thresh(mask, 
                            w = as.numeric(df_param$feature_op_threshold_w), 
                            h = as.numeric(df_param$feature_op_threshold_h), 
                            offset = as.numeric(df_param$feature_op_threshold_offset))
                        mask <- opening(mask, kern=makeBrush(as.numeric(df_param$feature_op_brush_size), 
                            shape=as.character(df_param$feature_op_brush_shape)))
                        seeder <- bwlabel(mask)
                        mask <- fillHull(thresh(mask, 
                            w = as.numeric(df_param$feature_op_brush_w), 
                            h = as.numeric(df_param$feature_op_brush_h), 
                            offset = as.numeric(df_param$feature_op_brush_offset)))
                        mask <- propagate(cell_stack_list[,,2,cell_idx], seeder, mask=mask)
                    
                    message(paste0("Number of features for ", cdat$FeatureID[cell_idx],":"))
                    max(mask)
                    
                    if(max(mask) != 0){

                                a1 <- as.data.frame(computeFeatures.shape(mask, 
                                                                             cell_stack_list[,,2, cell_idx]))
                                b1 <- as.data.frame(computeFeatures.moment(mask, 
                                                                            cell_stack_list[,,2, cell_idx]))
                                c1 <- as.data.frame(computeFeatures.haralick(mask, 
                                                                             cell_stack_list[,,2, cell_idx]))
                                d1 <- as.data.frame(computeFeatures.basic(mask, 
                                                                            cell_stack_list[,,2, cell_idx]))

                                # Stack each ERBs
                                edat[[cell_idx]] <- cbind.data.frame(a1, b1, c1, d1)
                                edat[[cell_idx]]$eFeatureID <- paste0(cdat$FeatureID[cell_idx], "_",
                                                                    1:max(mask))
                                edat[[cell_idx]]$FeatureID <- cdat$FeatureID[cell_idx]
                                artifacts = max(mask)

                            }else {

                                edat[[cell_idx]] <- NULL
                                artifacts = 0
                            }

                    cdat$nfeatures[cell_idx] <- artifacts

                } # For the for loop

                    # Store the ERB data
                    edat <- compact(edat)
                    if(length(edat) > 0){
                        edata <- do.call(rbind.data.frame, edat)
                        write.table(edata, paste(output, "/feature_matrix/feature-", 
                                             sample_id,".txt", 
                                            sep = ""), sep = "\t")
                        dat$edata <- edata

                    }else {
                        message("None features detected for the specified parameters")
#                         edata <- NULL
                    }

                

                    # Store the cell features
                    for(i in 1:nrow(cdat)){

                    # Store Cell features
                    imc <- cell_stack_list[,,,i]
                    writeImage(imc, paste("./images/segmented_images/stacked/cells/cell_", 
                                         cdat$FeatureID[i], ".png", 
                                          sep = ""), quality = 100)
                    }   

            }else {
                messsage("NONE Cells detected after rarefaction.")
            }

                write.table(cdat, 
                            paste(output, "/feature_matrix/cell-", 
                                        sample_id,".txt", 
                                        sep = ""), sep = "\t")
                dat$cdata <- cdat

        }
        
        # Return Object
        dat <- list(dat)
        names(dat) <- sample_id
        return(dat)
    } else{
        return(NULL)
    }

}, mc.cores = 27)
names(res) <- fls
res <- compact(res)
message(paste("Seggmentation Successfully executed!", sep = ""))
# Subset parameter 2 and compute features
glimpse(res)

# Create metadata
metadata <- data.frame(SampleID = as.character(fls))
metadata$SampleID <- str_replace(metadata$SampleID, 
                                 pattern = ".tif", replacement = "")
metadata$SampleID <- str_replace(metadata$SampleID, 
                                 pattern = as.character(live_stock), replacement = "")
metadata$SampleID <- str_replace(metadata$SampleID, 
                                 pattern = "/", replacement = "")

# Store procesed image binaries in Rdata format
res$metadata <- metadata

# Save image objects in an R container
save(list = "res", file = paste(data, "feature_image_data.Rdata", sep = ""))

# Write metadata in text format for plotting figures
write.table(metadata, paste(output, "metadata.txt", sep = ""), sep = "\t")

message(paste("---", "==Segmentation=="," Job ended at ", as.character(Sys.time()), "---", sep 
= ""))

# END OF SCRIPT
sessionInfo()
