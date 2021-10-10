#!/usr/bin/env Rscript
# Script for segmented cell wise analysis
# @Arpan Kumar Basak


rm(list = ls())

# Run this script for annotating the grayscale images only
# For colour donot use thresholding, instead treat the red channel as potential cellwall eroding the lines
message(paste("---", "==Computing Feature Statistics :: Image-wise analysis =="," Job Started at ", as.character(Sys.time()), "---", sep = ""))

# Script
rm(list = ls())
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
path <- as.character(getwd())
setwd(path)

# Load Requisites
pkgs <- c("tidyverse", "parallel", "vegan", "cowplot", "gridExtra", "patchwork")
lapply(pkgs, require, character.only = T)
options(warn = 1, mc.cores = 8)

# Measures not to fry R
options(set.cores = 8)
param <- "./script/parameters.R"
source(param)

figs <- paste0(figs, "cell_wise_analysis/")
stats <- paste0(stats, "cell_wise_analysis/")
output <- paste0(output, "cell_wise_analysis/")

dir_list <- c(figs, stats, output)

mclapply(dir_list, function(x){

  if(!dir.exists(x)){

  message("Directory created !!")

  dir.create(path = x, recursive = TRUE)

} else{

  message("Directory exists !!")
}
}, mc.cores = 4)

load(paste(data, "Obj_image_features.RData", sep = ""))


hmat <- as.matrix(imat_e[,grep(colnames(imat_e), pattern = "^h")])
smat <- as.matrix(imat_e[,grep(colnames(imat_e), pattern = "^s")])
mmat <- as.matrix(imat_e[,grep(colnames(imat_e), pattern = "^m")])
bmat <- as.matrix(imat_e[,grep(colnames(imat_e), pattern = "^b\\.")])

loc_metadata <- ifelse(!is.na(args[1]), as.character(args[1]), "./data/sample_data.txt")
k_clust <- ifelse(!is.na(args[2]), as.numeric(args[2]), 60)
k_optim <- args[3]

# Read metadata
meta <- read.delim2(loc_metadata, sep = "\t", header = T, as.is = T)

# Preparing tidy metadata
# meta$SampleID <- str_replace(meta$SampleID, "live_stock/", replacement = "")
# meta$SampleID <- str_replace(meta$SampleID, ".tif", replacement = "")
head(meta)

cname <- facts

# Make 
df <- cbind.data.frame(FeatureID = row.names(mmat), hmat, bmat, smat) %>%
    mutate(FeatureID = str_replace_all(FeatureID, "feature-", "")) %>%
    separate(FeatureID, into = c("SampleID", "x"), sep = "\\.", remove = FALSE, convert = FALSE) %>%
    cbind.data.frame(., meta[match(.$SampleID, meta$SampleID), which(colnames(meta) %in% facts)])

# Extract out morphology features eccentricity
levels(meta$genotype) <- as.character(genotype$name)
idx <- which(genotype$name %in% df$genotype)
geno.label <- geno.label[idx]
genotype <- genotype[idx,]

df$genotype <- factor(df$genotype, levels = genotype$name)
df$batch <- factor(df$batch, levels = batch$name)
meta$genotype <- factor(meta$genotype, levels = genotype$name)
df$SampleID <- as.factor(df$SampleID)
# df$batch <- as.factor(df$batch)


# Summarise
n_genotype <- df %>% group_by(genotype, SampleID) %>%
summarise(count = n()) %>% 
data.frame(.)

n_genotype_cells <- n_genotype %>% group_by(genotype) %>%
summarise(count = n()) %>% 
data.frame(.)

# mod <- aov(lm(formula = count ~ 0 + genotype:batch, data = n_genotype))
# fit <- broom::tidy(TukeyHSD(mod)) %>% na.omit(.) %>% mutate(significant = ifelse(adj.p.value <= 0.05, "*", "")) %>% arrange(desc(significant)) %>% data.frame(.)


n_genotype$n_images <- n_genotype_cells$count[match(n_genotype$genotype, n_genotype_cells$genotype)]

temp <- n_genotype %>% group_by(genotype) %>%
summarise(count = sum(count), n_images = sum(n_images))

write.table(temp, "./output/summary_experiment.txt", sep = "\t", quote = FALSE)

# Sample operations normalize matrix
df_samples <- cbind.data.frame(FeatureID = df$FeatureID, 
                            SampleID = df$SampleID, 
                            hmat, smat, bmat) %>%
            gather(key = "features", value = "vals", convert = FALSE, -FeatureID, -SampleID) %>%
            mutate(vals = as.numeric(vals)) %>%
            group_by(SampleID, features) %>%
            summarise(vals = mean(vals)) %>%
            spread(key = "features", value = "vals", convert = FALSE) %>%
            data.frame(.)
                            
row.names(df_samples) <- df_samples$SampleID
samples <- apply(df_samples[,-1], 2, function(x) (x - mean(x))/sd(x))

d <- 1 - cor(t(samples[,-1]))
mds_samples <- cmdscale(as.dist(d), eig = T, k = 3)
eigen <- round(100 * mds_samples$eig/sum(mds_samples$eig), 2)
hc <- hclust(as.dist(d), "average")
samples <- row.names(samples)[hc$order]
             
# Cluster features having the appropriate mean of the norphological features      
df_features <- cbind.data.frame(FeatureID = df$FeatureID, hmat, smat, bmat) %>%
            gather(key = "features", value = "vals", convert = FALSE, -FeatureID) %>%
            mutate(vals = as.numeric(vals)) %>%
            group_by(FeatureID, features) %>%
            summarise(vals = mean(vals)) %>%
            spread(key = "features", value = "vals", convert = FALSE, fill = 0) %>%
            data.frame(.)

# Argument
k <- k_clust

message("DONE..!!")

message("Initiating feature correlation...!!")

features <- apply(df_features[,-1], 2, as.numeric)
row.names(features) <- df_features$FeatureID
  
features <- apply(features, 2, function(x) (x - mean(x))/sd(x))

message("Computing optimum K value..!!")
n_centers=c(2, k_clust)
aic_bic_list <- mclapply(n_centers[1]:n_centers[2], function(x){

  set.seed(1)
  km_obj <- kmeans(features, centers=x, iter.max = 100)
  aic_bic <- k_optim(km_obj)

  #hc_clusters <- cutree(clust_obj, k=x)
 
  }

)
names(aic_bic_list) <- paste0("k_", as.character(n_centers[1]:n_centers[2]))
aic_bic_df <- do.call(rbind, aic_bic_list)
# kx <- last(which(aic_bic_df[,2]/sum(aic_bic_df[,2]) > 0.015))
cut_off = 0.0095
kx <- last(which(aic_bic_df[,1]/sum(aic_bic_df[,1]) > cut_off))
aic_bic_df <- as.data.frame(aic_bic_df) %>%
mutate(optim_AIC = AIC/sum(AIC), clusters = as.numeric(str_replace_all(row.names(.), "k_", "")))

# Plot the K optimum
ggplot(aic_bic_df, aes(x = clusters, y = optim_AIC)) +
geom_point(size = 1, colour = "black") +
geom_line(lty = "solid", lwd = 0.6) +
geom_vline(xintercept = kx, lty = "dashed", lwd = 0.3, colour = "darkred")+
geom_hline(yintercept = cut_off, lty = "dashed", lwd = 0.3, colour = "darkred")+
theme_bw() +
        labs(x = "Range of K based on number of genotypes",
             y = "AIC/sum(AIC)") +
       ggsave(filename = "./figures/optimised_kmeans.png", 
                              width = 4, height = 3, units = "in", 
              device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")

                  
message("Computing Correlation between features..!!")
d <- 1-cor(t(features))

# Kmeans clustering for 5 x 5 clusters
set.seed(1)
k=iflelse(!is.na(k_optim), kx, k_clust)
km <- kmeans(features, centers = k)
# spc <- kernlab::specc(features, centers = k)
hc <- hclust(as.dist(d), "average")
features_clusters <- row.names(features)[hc$order]
mds_cells <- cmdscale(as.dist(d), eig = T, k = 3)
eigen <- round(100 * mds_cells$eig/sum(mds_cells$eig), 2)


df_features_hmap <- cbind.data.frame(features, FeatureID = df_features[,1]) %>%
gather(key = "feature", value = "vals", convert = FALSE, -FeatureID) %>%
mutate(
                    FeatureID = factor(FeatureID, levels = features_clusters), 
                       feature = as.factor(feature), 
                       genotype = df$genotype[match(.$FeatureID, df$FeatureID)],
                       # batch = df$batch[match(.$FeatureID, df$cellID)],
                       vals = vals
                       ) %>%
                  # mutate(clusters = km$cluster[.$FeatureID]) %>%
data.frame(.)

# Make Profile facets
df_features_hmap$profile <- ""
df_features_hmap$profile[str_detect(as.character(df_features_hmap$feature), "^h")] <- "Haralick"
df_features_hmap$profile[str_detect(as.character(df_features_hmap$feature), "^b")] <- "Fluroscence Intensity"
df_features_hmap$profile[str_detect(as.character(df_features_hmap$feature), "^s")] <- "Spatial"
  

# Make MDS dataframe
mds_df <- data.frame(FeatureID = row.names(mds_cells$points), 
                       mds1 = mds_cells$points[,1], 
                       mds2 = mds_cells$points[,2], 
                       mds3 = mds_cells$points[,3]
                      ) %>%
cbind.data.frame(., df[match(.$FeatureID, df$FeatureID), which(colnames(df) %in% facts)])

mds_df$genotype <- factor(mds_df$genotype, levels = genotype$name[genotype$name %in% mds_df$genotype])
mds_df$batch <- factor(mds_df$batch, levels = batch$name[batch$name %in% mds_df$batch])

# Mean data frame
mds_df_mean <- mds_df %>%
group_by(genotype) %>%
summarise(mmds1 = mean(mds1), mmds2 = mean(mds2), mmds3 = mean(mds3)) %>%
data.frame(.)

(mds1_2 <- mds_df %>% 
     ggplot(aes(x= saturate(mds1), 
               y = saturate(mds2))) +
        geom_point(aes(colour = genotype), alpha = 0.4, size = 0.8) +
        geom_point(data = mds_df_mean, 
          aes(x = mmds1, y = mmds2, 
        colour = genotype), alpha = 0.8, size = 3) +
        geom_hline(yintercept = 0.0, alpha = 0.5, 
                colour = "darkgrey", 
                lty = "dashed") +
        geom_vline(xintercept = 0.0, alpha = 0.5, 
                colour = "darkgrey", 
                lty = "dashed") +
        geom_rug(sides = "tr", aes(colour = genotype), position = "identity", alpha = 0.5) +
        scale_colour_manual(values = as.character(genotype$`col.idx`), 
                                label = geno.label) +
        theme_AKB +
        theme() +
        labs(x = paste("MDS-1: ", eigen[1], " %", sep = ""),
             y = paste("MDS-2: ", eigen[2], " %", sep = ""),
             colour = "", 
             size = "")) +
       ggsave(filename = paste(figs, "MDS_feature_plots_12.png", sep = ""), 
                              width = 9, height = 7, 
              device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")

 (mds2_3 <- mds_df %>% 
            ggplot(aes(x= saturate(mds3), 
                  y = saturate(mds2))) +
           geom_point(aes(colour = genotype), alpha = 0.4, size = 0.8) +
           geom_point(data = mds_df_mean, 
            aes(x = mmds3, y = mmds2, colour = genotype), alpha = 0.8, size = 3) +
           geom_hline(yintercept = 0.0, alpha = 0.5, 
                   colour = "darkgrey", 
                   lty = "dashed") +
           geom_vline(xintercept = 0.0, alpha = 0.5, 
                   colour = "darkgrey", 
                   lty = "dashed") +
           geom_rug(sides = "tr", aes(colour = genotype), position = "identity", alpha = 0.5) +
           scale_colour_manual(values = as.character(genotype$`col.idx`), 
                                   label = geno.label) +
           theme_AKB +
           theme() +
           labs(x = paste("MDS-3: ", eigen[3], " %", sep = ""),
                y = paste("MDS-2: ", eigen[2], " %", sep = ""),
                colour = "", 
                size = "")) +
         ggsave(filename = paste(figs, "MDS_feature_plots_23.png", sep = ""), 
                                width = 9, height = 7, 
                device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")



bxp <- mclapply(colnames(mds_df)[str_detect(colnames(mds_df), "^mds")], function(i) {

  # i = "mds1"
  temp <- mds_df
  colnames(temp)[colnames(temp) == i] <- "comp"
  ax <- temp %>% ggplot(aes(x = genotype, y = comp)) +
                 geom_boxplot(fill = NA, 
                 aes(colour = genotype), 
                 alpha = 0.8, 
                 outlier.colour = NA,
                 size = 0.3
                 ) +
                 scale_colour_manual(
                  values = as.character(genotype$`col.idx`), 
                               label = geno.label) +
                 theme_void() +
                  labs(x = "", y = "", colour = "")

  return(ax)


 }, mc.cores = 8)
names(bxp) <- colnames(mds_df)[str_detect(colnames(mds_df), "^mds")]

# features_clusters <- row.names(features)
# idx <- paste0("erb-", df$FeatureID)
k_cluster <- km$cluster
names(k_cluster) <- str_replace(names(k_cluster), "\\.", "-")

message("DONE..!!")

message("Correlation between features..!!")
# Feature correlations
cc_features <- as.data.frame(d) %>% 
                 add_column(source = row.names(.), .before = 1) %>%
                 gather(key = "sink", value = "pcc", convert = FALSE, -source) %>%
                 filter(source != sink) %>%
                 mutate(
                      x_genotype = df$genotype[match(.$source, as.character(df$FeatureID))],
                      y_genotype = df$genotype[match(.$sink, as.character(df$FeatureID))],
                      source = factor(source, levels = features_clusters),
                      sink = factor(sink, levels = features_clusters),
                      pcc = 1 - pcc
                     ) %>%
               mutate(x_genotype = factor(x_genotype, levels = genotype$name),
                      y_genotype = factor(y_genotype, levels = genotype$name)) %>%
               data.frame(.)

cc_features <- cc_features %>% 
mutate(x_genotype = factor(x_genotype, levels = genotype$name), 
  y_genotype = factor(y_genotype, levels = genotype$name)
  # x_batch = as.factor(x_batch),
  # y_batch = as.factor(y_batch)
  ) 

# Number of correlative features within cells
# cc_mat <- cc_features %>% 
# group_by(x_genotype, source, y_genotype) %>%
# summarise(
#   similar = sum(sign(pcc) == 1), 
#   variant = sum(sign(pcc) == -1)) %>%
# ungroup() %>%
# mutate(similar = similar/(similar + variant),
#   variant = variant/(similar + variant)) %>%
# gather(key = "relation", value = "RA", convert = FALSE, similar, variant) %>%
# mutate(relation = as.factor(relation)) %>%
#   data.frame(.)

# levels(cc_mat$y_genotype) <- genotype$short[which(genotype$name %in% cc_mat$y_genotype)]

# # Comparison of proportion of cells being similar and dissimilar
# mod_var <- aov(lm(RA ~ 0 + y_genotype, cc_mat[cc_mat$relation == "variant",]))
# fit_var <- broom::tidy(TukeyHSD(mod_var)) %>%
# arrange(adj.p.value) %>%
# data.frame(.)

# mod_sim <- aov(lm(RA ~ 0 + y_genotype, cc_mat[cc_mat$relation == "similar",]))
# fit_sim <- broom::tidy(TukeyHSD(mod_sim)) %>%
# arrange(adj.p.value) %>%
# data.frame(.)

# Correlation matrix
# (hmap0 <- cc_mat %>% 
#   mutate(logra = ifelse(!is.finite(log10(RA)), NA, log10(RA))) %>%
#   ggplot(aes(x = y_genotype, y = source)) +
#                    geom_raster(aes(fill = logra)) +
#                    facet_grid(x_batch + x_genotype ~ relation, space = "fixed", scales = "free", switch = "both") +
#                    scale_fill_gradient(low = "black", high = "yellow", na.value = "black",
#                     breaks = c(-4, -3, -2, -1, 0), limit = c(-5, 0),
#                     labels = c("0.01", "0.1", "1", "10", "100")) +
#                    theme_AKB +
#                    theme(axis.text.x = element_blank(),
#                          axis.text.y = element_blank(), 
#                          strip.text.x = element_text(angle = 90), 
#                          strip.text.y = element_text(angle = 180), 
#                          legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
#                    labs(x = "", y = "", fill = "RA")) +
#                    ggsave(filename = paste(figs, "sample_corRA_heatmap.png", sep = ""), 
#                                 width = 9, height = 14, device = "png", dpi = 600, 
#                           limitsize = FALSE, bg = "transparent")

# Plot feature featre correlation
(hmap1 <- cc_features %>% 
 mutate(source = factor(source, levels = features_clusters),
   sink = factor(sink, levels = (features_clusters))
   ) %>%
 ggplot(aes(x = source, y = sink)) +
 geom_raster(aes(fill = pcc), alpha = 1) +
 facet_grid(y_genotype ~ x_genotype, space = "free_x", scales = "free", switch = "both") +
 scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 0.0, 
                                           mid = "white", na.value = "darkgrey") +
                      theme_AKB +
                      theme(panel.spacing = unit(0.1, "lines"),
                       axis.text.x = element_blank(),
                            axis.text.y = element_blank(), 
                            strip.text.y = element_text(angle = 180), 
                            strip.text.x = element_text(angle = 90), 
                            legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
                      labs(x = "", y = "", fill = "PCC")) +
                     ggsave(filename = paste(figs, "feature_cells_clustered_heatmap_genotype.png", 
                                             sep = ""), 
                            width = 9,
                            height = 14, 
                            device = "png", dpi = 600, 
                            limitsize = FALSE, bg = "transparent")

# Heatmap of features
(hmap2 <- df_features_hmap %>% 
 # filter(!is.nan(vals)) %>%
 mutate(FeatureID = factor(FeatureID, levels = features_clusters),
   profile = as.factor(profile)) %>%
 ggplot(aes(x = feature, y = FeatureID)) +
                  geom_raster(aes(fill = ifelse(!is.na(vals), saturate(vals), NA)), na.rm = FALSE) +
                  facet_grid(genotype ~ profile, space = "free_x", scales = "free", switch = "both") +
                  scale_fill_gradient2(low = "deeppink", high = "darkcyan", na.value = "darkgrey") +
                  theme_AKB +
                  theme(panel.spacing = unit(0.1, "lines"),
                       axis.text.x = element_text(size = 8),
                        axis.text.y = element_blank(), 
                        strip.text.x = element_text(angle = 90), 
                        strip.text.y = element_text(angle = 180), 
                        legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
                  labs(x = "", y = "", fill = "Z-score")) +
                 ggsave(filename = paste(figs, "/cell_feature_heatmap.png", sep = ""), 
                              width = 9, height = 14, device = "png", dpi = 600, 
                        limitsize = FALSE, bg = "transparent")


# (hmap3 <- df_features_hmap %>% 
#  # filter(!is.nan(vals)) %>%
#  mutate(FeatureID = factor(FeatureID, levels = features_clusters),
#    profile = as.factor(profile), cluster = factor(k_cluster[.$FeatureID], levels = 1:k),
#    vals = ifelse(!is.nan(vals), saturate(vals), NA)) %>%
#  ggplot(aes(x = feature, y = FeatureID)) +
#                   geom_raster(aes(fill = vals)) +
#                   facet_grid(cluster + genotype ~ profile, space = "free_x", scales = "free", switch = "both") +
#                   scale_fill_gradient2(low = "deeppink", high = "darkcyan", na.value = "darkgrey") +
#                   theme_AKB +
#                   theme(panel.spacing = unit(0.1, "lines"),
#                        axis.text.x = element_text(size = 8),
#                         axis.text.y = element_blank(), 
#                         strip.text.x = element_text(angle = 90), 
#                         strip.text.y = element_text(angle = 180), 
#                         legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
#                   labs(x = "", y = "", fill = "Z-score")) +
#                  ggsave(filename = paste(figs, "/cell_feature_clustered_heatmap.png", sep = ""), 
#                               width = 7, height = 14, device = "png", dpi = 600, 
#                         limitsize = FALSE, bg = "transparent")

alpha <- 0.01

# Make cluster summary

# Clustering genotypes devoid of features resembling nai1 genotype
clusters_genotype <- df %>% mutate(clusters = factor(k_cluster[.$FeatureID], levels = 1:k)) %>%
group_by(clusters, genotype) %>%
summarise(nrep = n()) %>% 
spread(key = genotype, value = nrep, fill = 0, convert = FALSE) %>%
filter(`nai1-gfph` == 0) %>%
mutate(cluster = paste0("cluster", clusters)) %>%
ungroup() %>%
select(cluster) %>%
unlist(., use.names = FALSE)

# Making proortions of cluster representatives
cluster_df <- df %>% mutate(clusters = factor(k_cluster[.$FeatureID], levels = 1:k)) %>%
group_by(clusters, genotype) %>%
summarise(nrep = n()) %>% 
spread(key = clusters, value = nrep, fill = 0, convert = FALSE) %>%
data.frame(.)

d <- apply(cluster_df[, str_detect(colnames(cluster_df), "^X")], 2, function(x) (x)/sum(x))
hc <- hclust(as.dist(1-cor(d)), "average")
sorted_cluster <- c(colnames(cluster_df)[str_detect(colnames(cluster_df), "^X")])[hc$order]
sorted_cluster <- str_replace_all(sorted_cluster, "X", "")

cluster_df <- cbind.data.frame(genotype = cluster_df$genotype, 
  apply(cluster_df[, str_detect(colnames(cluster_df), "^X")], 2, function(x) (x)/sum(x))
  ) %>% gather(key = "cluster", value = "nrep", convert = FALSE, -genotype) %>%
  mutate(clusters = factor(str_replace_all(cluster, "X", ""), levels = sorted_cluster)) %>%
  data.frame(.)

write.table(cluster_df, "./statistics/cell_cluster_table.txt", sep = "\t", quote = FALSE)

# Barplot for cluster representatives
(barp <- cluster_df %>%
  mutate(clusters = factor(as.character(clusters), 
    levels = sorted_cluster)) %>%
 ggplot(aes(x = clusters, y = nrep)) +
 geom_bar( 
           aes(fill = genotype), stat = "identity", alpha = 0.8, size = 0.3, position = "stack"
 ) +
     # geom_jitter(aes(colour = genotype), alpha = 0.1, size = 0.1, shape = 3) + 
     # facet_grid(batch ~ ., space = "fixed", scale = "fixed", switch = "x") + 
     scale_fill_manual(values = as.character(genotype$col.idx), na.value = "lightgrey",
                         label = geno.label) + 
     theme_AKB + 
     theme(axis.text.x = element_blank(),
           strip.text.x = element_text(angle = 90),
           axis.text.y = element_text(face = "bold", 
                                                      size = 14, 
                                                      angle = 30, 
                                                      hjust = 1, 
                                                      vjust = 0.5, color = "darkgrey"),
                           legend.position = "top", legend.text = element_text(size = 10)
          ) +
                     labs(x = "",
                          y = "",
                          colour = "",
                          shape = "",
                          fill = "") + 
                     coord_flip()) +
                    ggsave(filename = paste(figs, 
                                            "/feature_cell_kmeans_rep_barplot.png", sep = ""), 
                           width = 7, height = 14, 
                           device = "png", dpi = 600, 
                           limitsize = FALSE, bg = "transparent")

# Cluster representatives for each cells
cluster_df <- df %>% mutate(clusters = factor(k_cluster[.$FeatureID], levels = 1:k)) %>%
group_by(clusters, genotype, SampleID, FeatureID) %>%
summarise(nrep = n()) %>% 
spread(key = clusters, value = nrep, fill = 0, convert = FALSE) %>%
data.frame(.)

# Accomodation of features in clusters for distinct genotype and cells
cluster_df <- cbind.data.frame(
  # batch = cluster_df[,1],
  genotype = cluster_df[,1], 
  SampleID = cluster_df[,2], 
  cellID = cluster_df[,3],
  apply(cluster_df[,-c(1,2,3)], 2, function(x) (x - mean(x))/sd(x))
  ) %>% 
gather(key = "cluster", value = "vals", convert= FALSE, -genotype, -cellID, -SampleID) %>%
mutate(cluster = factor(str_replace_all(cluster, "X", ""), levels = 1:k),
  cellID = factor(cellID, levels = features_clusters), 
  SampleID = factor(SampleID, levels = samples),
  ) 

# (cluster_ra <- cluster_df %>%
#  ggplot(aes(x = cluster, y = cellID)) +
#                   geom_raster(aes(fill = saturate(vals))) +
#                   facet_grid(batch + genotype ~ ., space = "free_x", scales = "free", switch = "both") +
#                   scale_fill_gradient(low = "deeppink", high = "darkcyan", na.value = "darkgrey") +
#                   theme_AKB +
#                   theme(panel.spacing = unit(0.1, "lines"),
#                        axis.text.x = element_text(size = 8),
#                         axis.text.y = element_blank(), 
#                         strip.text.x = element_text(angle = 90), 
#                         strip.text.y = element_text(angle = 180), 
#                         legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
#                   labs(x = "", y = "", fill = "log2(Relative Proportions)")) +
#                  ggsave(filename = paste(figs, "/cell_clusters_heatmap.png", sep = ""), 
#                               width = 7, height = 14, device = "png", dpi = 600, 
#                         limitsize = FALSE, bg = "transparent")


(cluster_zscore <- cluster_df %>%
  mutate(cluster = factor(cluster, levels = sorted_cluster)) %>%
 group_by(cluster, genotype, SampleID) %>%
 summarise(vals = sum(vals)) %>%
 ggplot(aes(y = cluster, x = SampleID)) +
                  geom_raster(aes(fill = saturate(vals))) +
                  facet_grid(.~ genotype, space = "free_x", scales = "free", switch = "both") +
                  scale_fill_gradient(low = "black", high = "yellow", na.value = "black") +
                  theme_AKB +
                  theme(panel.spacing = unit(0.1, "lines"),
                        axis.text.y = element_text(size = 10),
                        axis.text.x = element_blank(), 
                        strip.text.x = element_text(angle = 90), 
                        # strip.text.y = element_text(angle = 180), 
                        legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
                  labs(x = "", y = "", fill = "Z-score")) +
                 ggsave(filename = paste(figs, "/sample_clusters_heatmap.png", sep = ""), 
                              width = 12, height = 14, device = "png", dpi = 600, 
                        limitsize = FALSE, bg = "transparent")


message("Spilling Feature-feature Plot..!!")


df_features_bp <- df_features %>% 
                gather(key = "feature", value = "vals", convert = FALSE, -FeatureID) %>%
                add_column(FeatureIDs = str_replace_all(.$FeatureID, "\\.", "-")) %>%
                  add_column(genotype = df$genotype[match(as.character(.$FeatureID), df$FeatureID)]) %>%
                  # add_column(batch = df$batch[match(as.character(.$FeatureID), df$cellID)]) %>%
                  add_column(SampleID = df$SampleID[match(as.character(.$FeatureID), df$FeatureID)]) %>%
                  add_column(k_cluster = k_cluster[match(.$FeatureIDs, names(k_cluster))]
                  ) %>%
                  add_column(cluster_rep = table(k_cluster)[match(.$k_cluster, names(table(k_cluster)))]
                    ) %>%
                  mutate(FeatureIDs = factor(FeatureIDs, levels = unique(features_clusters)),
                         feature = as.factor(feature), 
                         vals = as.numeric(vals), 
                         SampleID = factor(SampleID, samples),
                         k_cluster = factor(paste0("cluster", k_cluster), levels = paste0("cluster", sorted_cluster))
                         ) %>%
                  filter(k_cluster %in% clusters_genotype) %>% 
                  data.frame(.)


alpha <- 0.05

# Strip of cluster representatives
df_features_bp %>% 
group_by(k_cluster, cluster_rep) %>% 
summarise(n = n()) %>%
select(k_cluster, cluster_rep) %>%
mutate(
  labs = str_replace(as.character(k_cluster), "cluster", ""),
  cluster_rep = cluster_rep/sum(cluster_rep)
  ) %>%
ggplot(aes(x = "", y = k_cluster)) +
geom_raster(aes(fill = cluster_rep)) +
geom_text(aes(label = as.character(labs)), size = 10, fontface = "bold", show.legend = FALSE) +
scale_fill_gradient(low = "lightgrey", high = "darkgrey",
  guide = FALSE
  ) +
theme_void() +
labs(x = "", y = "", fill = "") +
ggsave(filename = paste(figs, "/cluster_density.png", sep = ""), 
                              width = 2, height = 14, device = "png", dpi = 600, 
                        limitsize = FALSE, bg = "transparent")

# Cluster determination
levels(df_features_bp$genotype) <- genotype$short
df_features_bp_denoise <- df_features_bp %>% filter(genotype != "nai1") # Removing the background features

res <- mclapply(unique(as.character(df_features_bp_denoise$feature)), function(x){


    temp <- NULL
    ret <- list()
    df_temp <- NULL
    fet <- NULL
    pw.comp <- NULL
    temp <- df_features_bp_denoise[df_features_bp_denoise$feature == x,]
    fet <- as.character(unique(temp$k_cluster))
    i = NULL
    # Something is wrong here
    ret <- lapply(1:length(fet), function(i){
                #     i <-  "cluster113"; x <- "s.perimeter"
        df_temp <- temp[temp$k_cluster == fet[i],]
        message(paste0("Cluster = ", fet[i], "; Feature = ", x, "; Number of Genotype = ", length(unique(df_temp$genotype))))
        if(length(unique(df_temp$genotype)) > 1){

          

          fit <- aov(glm(formula = vals ~ 0 + genotype, data = df_temp))
          pw.comp <- broom::tidy(TukeyHSD(fit)) %>% 
          # separate(comparison, into = c("x_genotype", "y_genotype"), sep = "-", convert = FALSE) %>%
          mutate(significant = ifelse(adj.p.value < alpha, "*", "NS")) %>% 
          arrange(significant) %>% 
          # mutate(cluster = x, feature = fet[i]) %>%
          mutate(feature = x, cluster = fet[i]) %>%
          data.frame(.)


        }else {

          
          message(paste0("Not all representatives fulfill the genotypes at cluster ", x))
        }
      df_temp <- NULL
      return(pw.comp)
  })
if(!is.null(ret))
  {
    ret <- do.call(rbind.data.frame, ret)
  }
return(ret)
}, mc.cores = 8)

stat_df <- do.call(rbind.data.frame, res) %>%
mutate(FDR = p.adjust(adj.p.value, method = "fdr"), feature = as.factor(feature)) %>%
separate(comparison, into = c("A", "B"), sep = "-", convert = FALSE) %>%
mutate(A = factor(A, levels = genotype$short),
  B = factor(B, levels = genotype$short),
  profile = factor(str_split(feature, "\\.", simplify = TRUE)[,1], 
    labels = c("Intensity", "Haralick", "Spatial"), levels = c("b", "h", "s")), 
  cluster = factor(cluster, levels = paste0("cluster", sorted_cluster))
)

# Heatmap for pairwise analysis
stat_df %>% ggplot(aes(x = A, y = cluster, fill = -log10(adj.p.value))) +
geom_raster(alpha = 1) +
  geom_tile(fill = NA, aes(colour = significant), width = 0.9, height = 0.9, size = 0.5) +
  facet_grid(. ~ profile + feature + B , space = "free", scales = "free", switch = "both") +
  scale_fill_gradient(na.value = "darkgrey", 
    low = "white", 
    high = "red") +
  scale_colour_manual(values = c(`*` = "black", `NS` = "transparent"),
    labels = c(`*` = "FDR < 0.01", `NS` = "")
    ) +
  theme_AKB +
        labs(x = "", 
             y = "", 
             colour = "", fill = "log10[adj.p.value]") +
        theme(legend.title = element_blank(), 
              axis.title = element_text(hjust = 0.5, vjust = 0.5),
              axis.line = element_blank(),
              axis.text.x = element_text(angle = 90, size = 12, hjust = 0.5, vjust = 0.5),
              strip.text.x = element_text(angle = 90, size = 12, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5)
              ) +
        ggsave(paste(figs, "/cluster_features_pvalue_hmap.png", sep = ""), 
               dpi = 600, device = "png", 
               height = 14, 
               width = 28, limitsize = F, bg = "transparent")


# Overall cell-wise boxplots
df_features_hmap %>%
    mutate(SampleID = factor(FeatureID, levels = features_clusters)) %>%
    ggplot(aes(x = genotype, y = vals)) +
    geom_point(aes(colour = genotype), position = position_jitterdodge(), alpha = 0.1, size = 0.5, shape = 3) + 
    facet_wrap(feature ~ ., scales = "free", as.table = TRUE, ncol = 8, nrow = 10) +
    geom_boxplot(fill = "white", 
                 aes(colour = genotype), 
                 alpha = 0.8, 
                 outlier.colour = NA,
                 size = 0.3
                 ) + 
    scale_colour_manual(values = as.character(genotype$col.idx), 
                        label = geno.label) + 
    theme_AKB + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold", 
                                                     size = 12, 
                                                     angle = 30, 
                                                     hjust = 1, 
                                                     vjust = 0.5, color = "darkgrey"),
                          legend.position = "top", legend.text = element_text(size = 10)
         ) + 
                    labs(x = "",
                         y = "",
                         colour = "",
                         shape = "",
                         fill = "") + 
                    ggsave(filename = paste(figs, 
                                            "/feature_summary_plot.png", sep = ""), 
                           width = 18, height = 20, units = "in",
                           device = "png", dpi = 600, 
                           limitsize = FALSE, bg = "transparent")


df_features_bp$profile <- ""
df_features_bp$profile[str_detect(as.character(df_features_bp$feature), "^h")] <- "Haralick"
df_features_bp$profile[str_detect(as.character(df_features_bp$feature), "^b")] <- "Fluroscence Intensity"
df_features_bp$profile[str_detect(as.character(df_features_bp$feature), "^s")] <- "Spatial"

# Subset the clustered images
df_features_bp <- df_features_bp %>% spread(key = "feature", value = "vals", fill = 0, convert = FALSE) %>%
data.frame(.)

imat.clustered <- imat[which(imat$FeatureID %in% df_features_bp$FeatureID),]
imat.clustered$k_cluster <- df_features_bp$k_cluster[match(imat.clustered$FeatureID, df_features_bp$FeatureID)]
imat.clustered$cell_path <- as.character(paste0("cell_", imat.clustered$FeatureID, ".png"))
imat.clustered$genotype <- df_features_bp$genotype[match(imat.clustered$FeatureID, df_features_bp$FeatureID)]

# Spill canvas
aa <- bxp$mds1 + theme(legend.position = "none") + coord_flip()
bb <- bxp$mds3 + theme(legend.position = "none") + coord_flip()
cc <- bxp$mds2 + theme(legend.position = "none")

# MDS
a1 <- cowplot::plot_grid(NULL, aa, 
  NULL, 
  bb, NULL, 
  nrow=1, ncol = 5, 
  rel_height = c(0.1, 0.1, 0.1, 0.1, 0.1),
  rel_widths = c(0.2, 0.8, 0.2, 0.8, 0.2),
  scale = c(0.5)
  )

a2 <- plot_grid(
  mds1_2 + theme(legend.position = "none"), 
  cc,
  mds2_3 + labs(y = "") + theme(legend.position = "none", 
    axis.text.y = element_blank(), 
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank()),
  nrow = 1, ncol = 3, rel_widths = c(2, 0.4, 2), 
  scale = c(1, 0.8, 1)
  )

mds_composite <- cowplot::plot_grid(
  a1,
  a2, 
  get_legend(mds1_2),
  nrow = 3, ncol = 1, "tblr",
  rel_widths = c(1, 0.1, 1), 
  rel_height = c(.1, 4, .1),
  scale = c(1, 1, 0.5),
  greedy = FALSE
)

ggsave(mds_composite, filename = paste0(figs, "mds_feature_summary.pdf"), 
                           width = 12, height = 10, 
                           device = "png", dpi = 600, 
                           limitsize = FALSE, bg = "white")



# Heatmap
a1 <- cowplot::plot_grid(
  hmap1 + theme_void() + theme(legend.position = "none", 
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 180)), NULL,
  hmap2 + theme_void() + theme(legend.position = "none", 
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 180)), NULL, 
  barp + theme(legend.position = "none", axis.text.y = element_text(size = 8, angle=0)),
  nrow = 1, ncol = 5, 
  axis = "tblr", 
  rel_widths = c(2, 0.1, 2, 0.1, 2), 
  rel_height = c(1)
)

a2 <- cowplot::plot_grid(
  get_legend(hmap1), 
  get_legend(hmap2), 
  nrow = 1, ncol = 2, 
  align = "hv", 
  rel_widths = c(.1, .1), 
  rel_height = 1
)

plot_obj <- cowplot::plot_grid(
  a2, a1, get_legend(barp),
  nrow = 3, ncol = 1, 
  align = "hv", axis = "l",
  scale = c(0.8, 1), 
  greedy = FALSE
)


# Publication quality
mds1_2 + ggsave(filename = paste(figs, "PQ_feature_MDS12.png", sep = ""), 
                              width = 5, height = 5, 
                              units = "in", device = "png", dpi = 600, 
                        limitsize = FALSE, bg = "transparent")

mds2_3 + ggsave(filename = paste(figs, "PQ_feature_MDS23.png", sep = ""), 
                              width = 5, height = 5, 
                              units = "in", device = "png", dpi = 600, 
                        limitsize = FALSE, bg = "transparent")


# Each boxplots seperately
mclapply(unique(as.character(df_features_hmap$feature)), function(i){

  message(paste0("Spilling ", i))

  stat <- c(mn = mean(df_features_hmap$vals[df_features_hmap$feature == i]),
            md = median(df_features_hmap$vals[df_features_hmap$feature == i]),
            q05 = quantile(df_features_hmap$vals[df_features_hmap$feature == i], .05),
            q95 = quantile(df_features_hmap$vals[df_features_hmap$feature == i], .95))

  
  df_features_hmap %>%
      filter(feature == i) %>%
      mutate(FeatureID = factor(FeatureID, levels = samples)) %>%
      ggplot(aes(x = genotype, y = vals)) +
      geom_hline(yintercept = stat, 
        alpha = 0.5, 
                lwd = 0.8,
               colour = "darkgrey", 
               lty = "dashed") +
      geom_point(aes(colour = genotype), position = position_jitterdodge(), alpha = 0.1, size = 0.5, shape = 3) + 
      # facet_grid(. ~ batch, scales = "free_x") +
      geom_boxplot(fill = NA, 
                   aes(colour = genotype), 
                   alpha = 0.8, 
                   outlier.colour = NA,
                   size = 0.5
                   ) + 
      scale_colour_manual(values = as.character(genotype$col.idx), 
                          label = geno.label) + 
      theme_AKB + 
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(face = "bold", 
                                                       size = 10, 
                                                       angle = 30, 
                                                       hjust = 1, 
                                                       vjust = 0.5, color = "darkgrey"),
                            legend.position = "top", legend.text = element_text(size = 10)
           ) + 
                      labs(x = "",
                           y = "",
                           colour = "",
                           shape = "",
                           fill = "") + 
                      ggsave(filename = paste(figs, 
                                              "/individual/PQ_feature_(", i,")_feature_plot.png", sep = ""), 
                             width = 6, height = 4,unit = "in",
                             device = "png", dpi = 600, 
                             limitsize = FALSE, bg = "transparent")

}, mc.cores = 24)
  

ggsave(plot_obj, filename = paste0(figs, "feature_profile_summary.pdf"), 
                           width = 20, height = 30, 
                           device = "png", dpi = 600, 
                           limitsize = FALSE, bg = "white")

# Remove the clusters corresponding to background
imat.clustered <- imat.clustered %>% filter(genotype != "nai1")

save(list = "imat.clustered", file = paste0(data, "Clustered_feature_image_data.RData"))

message("DONE!!")
sessionInfo()
