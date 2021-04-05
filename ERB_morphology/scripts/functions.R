# Function to make the stitching of images
contrast_project <- function(i, fls, id) {
    
    # TODO::Add a conditional statement to check the channel wether gray or rgb
    # TODO::Include Messages at each step of analysis for log and progress
    message(paste("Processing file ", id[i], ""))
    #x <- (fls[idx])[11]
    #i <- 11
    set.seed(1)
    ix <- grep(fls, pattern = id[i], value = TRUE)
    
    # Load Image - set for ERBs and Cells - Fetch a concat image for erb and cell
    img <- readImage(unlist(ix), type = "tiff")
    
    if(dim(img)[2] > 1024) {
        message(paste0("Image Dimension: ", dim(img)[1], " x ", dim(img)[2]))
        img <-img[,1:1024,,]
    }
    
    if(dim(img)[1] > 1024){
        message(paste0("Image Dimension: ", dim(img)[1], " x ", dim(img)[2]))
        erb <- img[1:1024,,2,]
        cell <- img[1026:2049,,1,]
        colorMode(img) <- Grayscale
    } else {
        message(paste0("Image Dimension: ", dim(img)[1], " x ", dim(img)[2]))
        erb <- img[,,2,]
        cell <- img[,,1,]
        colorMode(img) <- Grayscale
    }
    
    # Conduct MaxContrastProjection
    cell <- contrastProjection(imageStack = cell, 
                             w_x = 3,
                             w_y = 3, 
                             brushShape = "box", smoothing = 5)
    
    erb <- contrastProjection(imageStack = erb, 
                             w_x = 1,
                             w_y = 1, 
                             brushShape = "box", smoothing = 5)
    
    # Output Image
    img <- rgbImage(red = cell, green = erb)[,,-3]
    colorMode(img) <- Grayscale
    img  <- normalize(img)
         
    message(paste("\nKnitting image ", id[i], sep = ""))
    writeImage(x = img, 
               files =  paste(live_stock, id[i], ".tif", sep = ""), 
               quality = 100, 
               type = "tiff", 
               bits.per.sample = 8)
    
    message("DONE!!")
  }


# Make image segmentatio
segment_cells <- function(ind, fls, output) {
    
    #param <- "/klaster/work/abasak/git_repo_erb/ERB_ImageAnalysis/pilot_cotyledon/script/parameters.R"
    #require("EBImage")
    #require("tidyverse")
    #source(param)
    #imat <- data.frame()
    #ind <- 110
    dat <- list()
    out = output
    x <- fls[ind]
    set.seed(1)
    
    # Load Image - set for ERBs
    message(paste("Processing file ", x, ""))
    img <- readImage(x, type = "tiff")
    if(validObject(img)){

        #colorMode(img) <- Grayscale
        #img <- EBImage::normalize(img)
        #img <- img/255

        brush <- makeBrush(size=51, shape = "gaussian", sigma=1)
        erb <- filter2(img[,,2], filter = brush) # Correct function by Gaussian blur before
        # Morphology operation for ERBs and cell

        # Extract ERB features
        message("Processing ER-Body like features")
        erb <- erb * (erb > otsu(erb))
#         erb <- erb * (erb > quantile(erb, 0.96))
        erb <- thresh(erb, w = 10, h = 10, offset = 0.0001)
        erb <- opening(erb, kern=makeBrush(3, shape="disc"))
        eseed <- bwlabel(erb)
        #erb <- fillHull(thresh(temp, w = 20, h = 20, offset = 0.00005))
        emask <- propagate(img[,,2], eseed, mask=erb)

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
        cell <- filter2(img[,,1], filter = brush)
        cell <- cell * (cell < quantile(cell, 0.95))
        cell <- fillHull(thresh(cell, w = 20, h = 20, offset = 0.000001))
        cell <- opening(cell, kern=makeBrush(3,shape="disc"))
        cell.o <- opening(cell < quantile(cell, 0.99))
        cell[cell.o>cell] <- cell.o[cell.o>cell]
        cell[vac > cell] <- vac[vac > cell]
        cell[erb > cell] <- erb[erb > cell]
        cmask <- propagate(img[,,1], vmask, lambda=1.0e-4, mask=cell)

        # Segmentation
        seg_obj <- paintObjects(vac, img[,,1], thick = T, col = "#ffff31")
        seg_obj <- paintObjects(cell, seg_obj, thick = F, col = "#ffff10")
        seg_obj <- paintObjects(erb, seg_obj, thick = F, col = "#ffff00")

        merge_obj <- rgbImage(red = img[,,1] * cell, 
                             green = img[,,2] * erb, 
                             blue = img[,,1] * vac)
        
        # Merge and paint the cell-like features and stack them
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
            id_fpc <- which(cdat[,1] > 10000) # Filter garbage or rarefy
            cdat$cell_feature <- 0
            cdat$cell_feature[id_fpc] <- 1
            if(length(id_fpc) != 0){
                
                cdat <- cdat %>% filter(cell_feature == 1) %>% add_column(SampleID = sample_id, .before = 1)
                cdat <- cdat %>% add_column(FeatureID = paste(sample_id, row.names(cdat), sep = "_"),
                                            .before = 1)
                cell_stack_list <- cell_stack_list[,,,id_fpc]
                message("Filtered potenial cells:")
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
                cdat$nerbs <- 0
                edat <- list()
                for(cell_idx in 1:nrow(cdat)){

                        # Segment ERbodies
        #                 cell_idx <- 13
                        mask <- cell_stack_list[,,2, cell_idx]
                        colorMode(mask) <- Grayscale
                        mask <- mask * (mask > otsu(mask))
#                         mask <- mask * (mask > quantile(mask, 0.9))
                        mask <- thresh(mask, w = 10, h = 10, offset = 0.0001)
                        mask <- opening(mask, kern=makeBrush(3, shape="disc"))
                        seeder <- bwlabel(mask)
                        mask <- fillHull(thresh(mask, w = 20, h = 20, offset = 0.005))
                        mask <- propagate(cell_stack_list[,,2,cell_idx], seeder, mask=mask)
                    
                    message(paste0("Number of ERBs for ", cdat$FeatureID[cell_idx],":"))
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

                    cdat$nerbs[cell_idx] <- artifacts

                } # For the for loop

                    # Store the ERB data
                    edat <- compact(edat)
                    if(length(edat) > 0){
                        edata <- do.call(rbind.data.frame, edat)
                        write.table(edata, paste(output, "/feature_matrix/erb-", 
                                             sample_id,".txt", 
                                            sep = ""), sep = "\t")
                        dat$edata <- edata

                    }else {
                        message("None ERB features")
#                         edata <- NULL
                    }

                

                    # Store the cell features
                    for(i in 1:nrow(cdat)){

                    # Store Cell features
                    imc <- cell_stack_list[,,,i]
                    writeImage(imc, paste("./images/segmented_images/stacked/cell_features/cell_", 
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
}

# 