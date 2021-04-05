#!/usr/bin/env Rscript
# Script for clustered cell feature analysis
# @Arpan Kumar Basak

# Script
rm(list = ls())
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
path <- as.character(getwd())
setwd(path)

# Run this script for annotating the grayscale images only
# For colour donot use thresholding, instead treat the red channel as potential cellwall eroding the lines
message(paste("---", "==Computing Feature Statistics=="," Job Started at ", as.character(Sys.time()), "---", sep = ""))

# Load Requisites
pkgs <- c("EBImage", "tidyverse", "parallel", "vegan", "mda", "klaR", "cowplot")
lapply(pkgs, require, character.only = T)
options(warn = 1, mc.cores = 8)


# Measures not to fry R
options(set.cores = 8)
param <- "./script/parameters.R"
source(param)

figs <- paste0(figs, "clustered_analysis/")
stats <- paste0(stats, "clustered_analysis/")
output <- paste0(output, "clustered_analysis/")

dir_list <- c(figs, stats, output)

mclapply(dir_list, function(x){

  if(!dir.exists(x)){

  message("Directory created !!")

  dir.create(path = x, recursive = TRUE)

} else{

  message("Directory exists !!")
}
}, mc.cores = 4)

cell_path <- paste0("./images/segmented_images/stacked/cells/")
feature_path <- paste0("./images/segmented_images/stacked/features/")

load(paste0(data, "Feature_Clustered_images_data.RData"))

# Isolate the headers
hmat <- as.matrix(imat.clustered[,grep(colnames(imat.clustered), pattern = "^h")])
smat <- as.matrix(imat.clustered[,grep(colnames(imat.clustered), pattern = "^s")])
mmat <- as.matrix(imat.clustered[,grep(colnames(imat.clustered), pattern = "^m")])
bmat <- as.matrix(imat.clustered[,grep(colnames(imat.clustered), pattern = "^b\\.")])


# Prerequisites
df_features <- cbind.data.frame(
  FeatureID = imat.clustered$eFeatureID, 
  k_cluster = imat.clustered$k_cluster,
                            hmat, smat, bmat) %>%
            gather(key = "features", value = "vals", convert = FALSE, -FeatureID, -k_cluster) %>%
            mutate(vals = as.numeric(vals)) %>%
            group_by(FeatureID, k_cluster, features) %>%
            summarise(vals = mean(vals)) %>%
            spread(key = "features", value = "vals", convert = FALSE, fill = 0) %>%
            data.frame(.)

# Make a metadata table to avoid this step
row.names(df_features) <- df_features$FeatureID
feature <- apply(df_features[,-c(1, 2)], 2, function(x) (x - mean(x))/sd(x))
df <- df_features %>% dplyr::select(FeatureID, k_cluster) %>% data.frame(.)
row.names(df) <- df$FeatureID



df$Genotype <- df$genotype
idx <- which(as.character(genotype$name) %in% df$Genotype)
df$Genotype <- factor(df$Genotype, levels = as.character(genotype$name[idx]))
df$date <- as.factor(str_split(df$FeatureID, pattern = "_", simplify = TRUE)[,1])
df$k_cluster <- as.factor(df$k_cluster)
df$batch <- as.factor(df$batch)
df$group <- as.factor(paste(df$Genotype, df$batch, sep = "_"))

# Make the data table as html file and include the clustered cells

# Define model formula
f <- formula(as.dist(d) ~ group + Condition(k_cluster))

d <- 1 - cor(t(feature))
hc <- hclust(as.dist(d), "average")
feature_clusters <- row.names(feature)[hc$order]
class_weights <-  100/table(df$Genotype)

# Conduct a CCA
if(file.exists("./data/cca_data.Rdata")){
  load("./data/cca_data.Rdata")
  } else{
   
   cca_obj <- vegan::capscale(formula = f, df, 
    distance = "pearson", sqrt.dist = T, add = F)
  
  # Conduct PERMANOVA
  set.seed(1)
  p.nova <- vegan::anova.cca(cca_obj, by = "term", parallel = 8)

    save(list = c("cca_obj", "p.nova"), file = "./data/cca_data.Rdata") 
  }

# Extract variables from stat objects
eig <- format(100 * (cca_obj$CCA$eig/sum(cca_obj$CCA$eig)), digits = 4)
pval <- p.nova$`Pr(>F)`[1]
chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                       row.names = c("total", "constrianed", "unconstrained"))
write.table(variable, "./statistics/cca_results.txt", sep = "\t", quote = F)

# Format P-Value for representation
variance <- format(100 * variable$proportion[2], digits = 4)
ti <- paste("Genotype = ", variance, "% p = ", pval, 
             ifelse(pval < 0.05, "*", " (NS)"), sep = "")

# Conduct MDA
mda_obj <- mda::mda(formula = formula(df$Genotype ~ as.matrix(feature)), 
  subclasses = length(unique(df$Genotype))
)

temp <- as.data.frame(confusion(mda_obj)) %>% 
spread(key = true, value = Freq, convert = FALSE)
mat <- apply(temp[,-1], 1, function(x) (x)/sum(x))
colnames(mat) <- colnames(temp)[-1]

# Plot Confusion Matrix
(cmap1 <- cbind.data.frame(pred = row.names(mat), mat) %>%
   gather(key = "true", value = "prop", convert = FALSE, -pred) %>%
   mutate(true = str_replace_all(true, "\\.", "-")) %>%
   mutate(pred = factor(pred, levels = as.character(genotype$name)),
     true = factor(true, levels = as.character(genotype$name)),
     prop = (100*prop)
     ) %>%
   ggplot(aes(x = pred, y = true, fill = prop)) +
   geom_raster(alpha = 1) +
   scale_fill_gradient(low = "white", high = "red", 
     breaks = c(0,25, 50, 75, 100), limits = c(0, 100),
     labels = paste0(seq(0,100,25), "%")) +
   theme_AKB +
         labs(x = "True Class", 
              y = "Predicted Class", 
              colour = "") +
         theme(legend.title = element_blank(), 
               axis.title = element_text(hjust = 0.5, vjust = 0.5),
               axis.line = element_blank(),
               axis.text.x = element_text(angle = 90, size = 12, hjust = 0.5, vjust = 0.5),
               axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5)
               )) +
        ggsave(paste(figs, "/feature_mda_cm.png", sep = ""), 
               dpi = 600, device = "png", 
               height = 7, 
               width = 7, limitsize = F, bg = "transparent")

# Conduct FDA
fda_obj <- mda::fda(formula = formula(df$Genotype ~ as.matrix(feature)), 
  weights = class_weights
)

temp <- as.data.frame(confusion(fda_obj)) %>% 
spread(key = true, value = Freq, convert = FALSE)
mat <- apply(temp[,-1], 1, function(x) (x)/sum(x))
colnames(mat) <- colnames(temp)[-1]

# Plot Confusion Matrix
(cmap2 <- cbind.data.frame(pred = row.names(mat), mat) %>%
  gather(key = "true", value = "prop", convert = FALSE, -pred) %>%
  mutate(true = str_replace_all(true, "\\.", "-")) %>%
  mutate(pred = factor(pred, levels = as.character(genotype$name)),
    true = factor(true, levels = as.character(genotype$name)),
    prop = (100*prop)
    ) %>%
  ggplot(aes(x = pred, y = true, fill = prop)) +
  geom_raster(alpha = 1) +
  scale_fill_gradient(low = "white", high = "red", 
    breaks = c(0,25, 50, 75, 100), limits = c(0, 100),
    labels = paste0(seq(0,100,25), "%")) +
  theme_AKB +
        labs(x = "True Class", 
             y = "Predicted Class", 
             colour = "") +
        theme(legend.title = element_blank(), 
              axis.title = element_text(hjust = 0.5, vjust = 0.5),
              axis.line = element_blank(),
              axis.text.x = element_text(angle = 90, size = 12, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5)
              ) )+
        ggsave(paste(figs, "/feature_fda_cm.png", sep = ""), 
               dpi = 600, device = "png", 
               height = 7, 
               width = 7, limitsize = F, bg = "transparent")

# Plot CCA -- find the location of the factors
marker_x <- sapply(unique(df$Genotype), function(x) 3 * median(as.numeric(cca_obj$CCA$wa[row.names(df),1])[df$Genotype == x]))
marker_y <- sapply(unique(df$Genotype), function(x) 3 * median(as.numeric(cca_obj$CCA$wa[row.names(df),2])[df$Genotype == x]))
idx <- which(as.character(genotype$name) %in% df$Genotype)
names(marker_x) <- genotype$name[idx]
names(marker_y) <- genotype$name[idx]


(cca_plot <- df %>% 
     dplyr::mutate(CCA1 = as.numeric(cca_obj$CCA$wa[row.names(df),1]),
                   CCA2 = as.numeric(cca_obj$CCA$wa[row.names(df),2])) %>% 
     ggplot2::ggplot(aes(x= (CCA1), y = (CCA2))) +
     ggtitle(ti) + 
     geom_point(alpha = 0.4, 
                size = 0.5, 
                aes(colour = Genotype)) +
     geom_hline(yintercept = 0.0, linetype = "dashed", colour = "darkgrey", alpha = 0.6, size = 0.8) +
     geom_vline(xintercept = 0.0, linetype = "dashed", colour = "darkgrey", alpha = 0.6, size = 0.8) +
         # geom_segment(aes(x = 0, y = 0, xend = marker_x[1], yend = marker_y[1]), color =  as.character(genotype$col.idx)[idx][1],
         #   arrow = arrow(length = unit(1, "cm"))) +
         # geom_segment(aes(x = 0, y = 0, xend = marker_x[2], yend = marker_y[2]), color =  as.character(genotype$col.idx)[idx][2],
         #   arrow = arrow(length = unit(1, "cm"))) +
         # geom_segment(aes(x = 0, y = 0, xend = marker_x[3], yend = marker_y[3]), color =  as.character(genotype$col.idx)[idx][3],
         #   arrow = arrow(length = unit(1, "cm"))) +
         # geom_segment(aes(x = 0, y = 0, xend = marker_x[4], yend = marker_y[4]), color =  as.character(genotype$col.idx)[idx][4],
         #   arrow = arrow(length = unit(1, "cm"))) +
         # geom_segment(aes(x = 0, y = 0, xend = marker_x[5], yend = marker_y[5]), color =  as.character(genotype$col.idx)[idx][5],
         #   arrow = arrow(length = unit(1, "cm"))) +
         # geom_segment(aes(x = 0, y = 0, xend = marker_x[6], yend = marker_y[6]), color =  as.character(genotype$col.idx)[idx][6],
         #   arrow = arrow(length = unit(1, "cm"))) +
         # geom_segment(aes(x = 0, y = 0, xend = marker_x[7], yend = marker_y[7]), color =  as.character(genotype$col.idx)[idx][7],
         #   arrow = arrow(length = unit(1, "cm"))) +
     stat_ellipse(type = "norm", alpha = 0.6, linetype = "dashed",
                  level = 0.95, aes(group = Genotype, 
                                    colour = Genotype)) +
     scale_color_manual(values = as.character(genotype$col.idx)[idx], labels = geno.label[idx]) + 
     # scale_linetype_manual(values = genotype$lty.idx, guide = FALSE) +
     theme_AKB +
     theme(legend.title = element_blank(), 
       axis.title = element_text(vjust = 0.5, hjust = 0.5, size = 14),
       axis.text = element_text(size = 12, vjust = 0.5, hjust = 0.5),
               panel.border = element_rect(size = 0.8, fill = "transparent"), 
               axis.line = element_blank()) +
     labs(x = paste("CCA 1: ",eig[1],"%", sep = ""), 
          y = paste("CCA 2: ",eig[2],"%", sep = ""), 
          colour = "", shape = "", lty = "")) +
    ggsave(filename = paste(figs,"/feature_cca_analysis.png", sep = ""), 
           dpi = 600, width = 10, height = 10, device = "png", 
           limitsize = FALSE, bg = "transparent")

# Based on the 3rd argument tSNE and spectral clusterig of organelle features

# Conduct tSNE clustering here for features
if(file.exists("./data/tsne_clusetered.Rdata")){
  load("./data/tsne_clusetered.Rdata")
  } else{
    
    set.seed(1)
    tsne_obj <- Rtsne::Rtsne(as.dist(d)^2, 
    dims = 2, 
    initial_dims = 60,
    perplexity = 30, 
    theta = 0.05, 
    check_duplicates = TRUE,
    partial_pca = TRUE, 
    max_iter = 1000,
    vfeatureose = getOption("vfeatureose", FALSE), 
    is_distance = TRUE,
    Y_init = NULL, 
    normalize = TRUE, 
    momentum = 0.05, 
    final_momentum = 0.8, 
    eta = 200,
    exaggeration_factor = 10, 
    num_threads = 24
  )

    save(list = "tsne_obj", file = "./data/tsne_clusetered.Rdata") 
}

# idx <- which(as.character(genotype$name) %in% df$Genotype)
# marker_x <- sapply(unique(df$Genotype), function(x) 3 * median(as.numeric(tsne_obj$Y[,1])[df$Genotype == x]))
# marker_y <- sapply(unique(df$Genotype), function(x) 3 * median(as.numeric(tsne_obj$Y[,2])[df$Genotype == x]))
# names(marker_x) <- genotype$name[idx]
# names(marker_y) <- genotype$name[idx]

# # =======
# median_arrows <- function(marker_x, marker_y) {
  
#   idx <- genotype$name == names(marker_x)
#   obj <- NULL
#   for(ix in 1:length(marker_x)){

#     id <- idx[ix]
#     if(is.null(obj)){
#       obj <- ggplot2::geom_segment(aes(x = 0, y = 0, xend = marker_x[ix], yend = marker_y[ix]), color =  as.character(genotype$col.idx)[idx][1],
#             arrow = arrow(length = unit(1, "cm"))) 
#     } else{
#       obj <- obj + ggplot2::geom_segment(aes(x = 0, y = 0, xend = marker_x[ix], yend = marker_y[ix]), color =  as.character(genotype$col.idx)[idx][1],
#             arrow = arrow(length = unit(1, "cm")))
#     }
    
#   }
#  return(obj)   
# }
# ===========


# ======

# Conduct SEPCC clustering
# if(file.exists("./data/specs_clustered.Rdata")){
#   load("./data/specs_clusetered.Rdata")
#   } else{
    
#     set.seed(1)
#     spec_obj <- kernlab::specc(x = d[1:10000,1:10000],
#     centers = 8,
#     kernel = "rbfdot", 
#     kpar = "automatic", 
#     iterations = 100, 
#     mod.sample = 0.75
#   )

#     save(list = "spec_obj", file = "./data/specs_clusetered.Rdata") 
# }

# ======

# df %>% 
#     dplyr::mutate(tSNE1 = as.numeric(tsne_obj$Y[,1]),
#                   tSNE2 = as.numeric(tsne_obj$Y[,2])) %>% 
#     ggplot(aes(x = (tSNE1), y = (tSNE2))) +
#         geom_point(alpha = 0.4, aes(colour = Genotype), size = 1) +
#         geom_hline(yintercept = 0.0, linetype = "dashed", colour = "darkgrey", alpha = 0.6, size = 0.8) +
#         geom_vline(xintercept = 0.0, linetype = "dashed", colour = "darkgrey", alpha = 0.6, size = 0.8) +
#         scale_color_manual(values = as.character(genotype$col.idx)[idx], labels = geno.label[idx]) + 
#         scale_fill_manual(values = as.character(genotype$col.idx)[idx], labels = geno.label[idx]) + 
#         stat_ellipse(type = "euclid", alpha = 0.6, linetype = "dashed",
#                  level = 0.95, aes(group = Genotype, 
#                                    colour = Genotype)) +
#         theme_AKB +
#         labs(x = paste("tSNE 1: "),
#              y = paste("tSNE 2: "),
#              colour = "",
#              shape = "") +
#         theme(legend.title = element_blank(), 
#           axis.title = element_text(hjust = 0.5, vjust = 0.5),
#               panel.border = element_rect(size = 0.8, fill = "transparent"), 
#               axis.line = element_blank()) +
#         ggsave(paste(figs, "/feature_tSNE_Plot_1.png", sep = ""), 
#                dpi = 600, device = "png", 
#                height = 10, 
#                width = 9, limitsize = F, bg = "transparent")


# (tsne_plot <- df %>% 
#      dplyr::mutate(tSNE1 = as.numeric(tsne_obj$Y[,1]),
#                    tSNE2 = as.numeric(tsne_obj$Y[,2])) %>% 
#      ggplot(aes(x = (tSNE1), y = (tSNE2))) +
#          geom_point(alpha = 0.4, aes(colour = Genotype), size = 1) +
#          geom_hline(yintercept = 0.0, linetype = "dashed", colour = "darkgrey", alpha = 0.6, size = 0.8) +
#          geom_vline(xintercept = 0.0, linetype = "dashed", colour = "darkgrey", alpha = 0.6, size = 0.8) +
#          scale_color_manual(values = as.character(genotype$col.idx)[idx], labels = geno.label[idx]) + 
#          scale_fill_manual(values = as.character(genotype$col.idx)[idx], labels = geno.label[idx]) + 
#          facet_grid(.~batch, space = "free", scale = "free", switch = "both") + 
#          geom_segment(aes(x = 0, y = 0, xend = marker_x[1], yend = marker_y[1]), color =  as.character(genotype$col.idx)[idx][1],
#            arrow = arrow(length = unit(1, "cm"))) +
#          geom_segment(aes(x = 0, y = 0, xend = marker_x[2], yend = marker_y[2]), color =  as.character(genotype$col.idx)[idx][2],
#            arrow = arrow(length = unit(1, "cm"))) +
#          geom_segment(aes(x = 0, y = 0, xend = marker_x[3], yend = marker_y[3]), color =  as.character(genotype$col.idx)[idx][3],
#            arrow = arrow(length = unit(1, "cm"))) +
#          geom_segment(aes(x = 0, y = 0, xend = marker_x[4], yend = marker_y[4]), color =  as.character(genotype$col.idx)[idx][4],
#            arrow = arrow(length = unit(1, "cm"))) +
#          geom_segment(aes(x = 0, y = 0, xend = marker_x[5], yend = marker_y[5]), color =  as.character(genotype$col.idx)[idx][5],
#            arrow = arrow(length = unit(1, "cm"))) +
#          geom_segment(aes(x = 0, y = 0, xend = marker_x[6], yend = marker_y[6]), color =  as.character(genotype$col.idx)[idx][6],
#            arrow = arrow(length = unit(1, "cm"))) +
#          geom_segment(aes(x = 0, y = 0, xend = marker_x[7], yend = marker_y[7]), color =  as.character(genotype$col.idx)[idx][7],
#            arrow = arrow(length = unit(1, "cm"))) +
#          theme_AKB +
#          labs(x = paste("tSNE 1: "),
#               y = paste("tSNE 2: "),
#               colour = "",
#               shape = "") +
#          theme(legend.title = element_blank(), 
#            axis.title = element_text(hjust = 0.5, vjust = 0.5),
#                panel.border = element_rect(size = 0.8, fill = "transparent"), 
#                axis.line = element_blank())) +
#         ggsave(paste(figs, "/feature_tSNE_Plot_3.png", sep = ""), 
#                dpi = 600, device = "png", 
#                height = 7, 
#                width = 14, limitsize = F, bg = "transparent")

# (tsne_plot + ggsave(paste(figs, "/PQ_feature_tSNE_Plot_3.png", sep = ""), 
#                dpi = 600, device = "png", units = "in",
#                height = 4, 
#                width = 8, limitsize = F, bg = "transparent"))


# # Spill canvas
# a1 <- cowplot::plot_grid(cmap1 + theme(legend.position = "none"), 
#   cmap2 + theme(legend.position = "none"), 
#   nrow=1, ncol = 2, 
#   rel_height = c(1),
#   rel_widths = c(1,1),
#   scale = c(0.8)
#   )

# a2 <- cowplot::plot_grid(get_legend(cmap1), 
#   a1,
#   cca_plot + theme(legend.position = "none"), 
#   nrow=3, ncol = 1, 
#   rel_height = c(0.1, 0.3, 1),
#   rel_widths = c(0.4, 1, 1),
#   scale = c(1)
#   )

# composite_features <- cowplot::plot_grid(
#   a2,
#   tsne_plot + theme(legend.position = "none"),
#   nrow = 1, ncol = 2, "tblr",
#   rel_widths = c(1, 1), 
#   rel_height = c(1),
#   scale = c(1, 1),
#   greedy = FALSE
# )

# ggsave(composite_features, 
#   filename = paste(figs, "/profile_feature_summary.png", sep = ""), 
#                dpi = 600, device = "png", 
#                height = 10, 
#                width = 20, limitsize = FALSE, bg = "transparent")


# composite_features <- cowplot::plot_grid(
#   a2,
#   specc_plot + theme(legend.position = "none"),
#   nrow = 1, ncol = 2, "tblr",
#   rel_widths = c(1, 1), 
#   rel_height = c(1),
#   scale = c(1, 1),
#   greedy = FALSE
# )

# ggsave(composite_features, 
#   filename = paste(figs, "/profile_feature_summary.png", sep = ""), 
#                dpi = 600, device = "png", 
#                height = 10, 
#                width = 20, limitsize = FALSE, bg = "transparent")

# Plot the results from Spectral cllustering

# Stats

# Plot Canvas


message("DONE!!")
sessionInfo()