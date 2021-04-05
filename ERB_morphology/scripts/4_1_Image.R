#!/usr/bin/env Rscript
# Script for sample wise analysis
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

figs <- paste0(figs, "sample_wise_analysis/")
stats <- paste0(stats, "sample_wise_analysis/")
output <- paste0(output, "sample_wise_analysis/")

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

# Read metadata
meta <- read.delim2(loc_metadata, sep = "\t", header = T, as.is = T)

# Preparing tidy metadata
meta$SampleID <- str_replace(meta$SampleID, "live_stock/", replacement = "")
meta$SampleID <- str_replace(meta$SampleID, ".tif", replacement = "")
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
sum_cells <- df %>% group_by(genotype) %>%
summarise(n = n()) %>%
data.frame(., stringsAsFactors = FALSE)

sum_images <- df %>% group_by(genotype, SampleID) %>%
summarise(n = n()) %>%
ungroup(.) %>%
group_by(genotype) %>%
summarise(n = n()) %>%
data.frame(., stringsAsFactors = FALSE)

write.table((cbind.data.frame(sum_cells, n_img = sum_images$n) %>% mutate(ratio = n/n_img)), 
  "./output/summary_all.txt", 
  quote = FALSE, row.names = FALSE, sep = "\t")


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
samples <- apply(df_samples[,-1], 2, function(x) (x-mean(x))/sd(x))

d <- 1 - cor(t(samples[,-1]))
mds_samples <- cmdscale(as.dist(d), eig = T, k = 3)
eigen <- round(100 * mds_samples$eig/sum(mds_samples$eig), 2)
hc <- hclust(as.dist(d), "average")
samples <- row.names(samples)[hc$order]

# Make MDS dataframe
mds_df <- data.frame(
                      SampleID = row.names(mds_samples$points), 
                       mds1 = mds_samples$points[,1], 
                       mds2 = mds_samples$points[,2], 
                       mds3 = mds_samples$points[,3]
                      ) %>%
cbind.data.frame(., meta[match(.$SampleID, meta$SampleID), which(colnames(meta) %in% facts)])
                 

# Sample correlations
cc_samples <- as.data.frame(d) %>% 
                 add_column(source = row.names(.), .before = 1) %>%
                 gather(key = "sink", value = "pcc", convert = FALSE, -source) %>%
                 filter(source != sink) %>%
                 mutate(x_genotype = meta$genotype[match(source, as.character(meta$SampleID))],
                        y_genotype = meta$genotype[match(sink, as.character(meta$SampleID))],
                        x_batch = meta$batch[match(source, as.character(meta$SampleID))],
                        y_batch = meta$batch[match(sink, as.character(meta$SampleID))],
                        source = factor(source, levels = samples),
                        sink = factor(sink, levels = samples),
                        pcc = 1 - pcc
                       ) %>%
                 mutate(x_genotype = factor(x_genotype, levels = genotype$name),
                        y_genotype = factor(y_genotype, levels = genotype$name),
                        x_batch = factor(x_batch, levels = batch$name[batch$name %in% .$x_batch]),
                        y_batch = factor(y_batch, levels = batch$name[batch$name %in% .$y_batch])) %>%
                 data.frame(.)

# Number of correlative features within cells
# cc_mat <- cc_samples %>% 
# # filter(as.character(x_batch) == as.character(y_batch)) %>%
# group_by(x_genotype, source, y_genotype) %>%
# summarise(similar = sum(sign(pcc) == 1), 
#   variant = sum(sign(pcc) == -1)) %>%
# ungroup() %>%
# mutate(similar = similar/(similar + variant),
#   variant = variant/(similar + variant)) %>%
# gather(key = "relation", value = "RA", convert = FALSE, similar, variant) %>%
# mutate(relation = as.factor(relation)) %>%
#   data.frame(.)

# levels(cc_mat$x_genotype) <- genotype$short[which(genotype$name %in% cc_mat$x_genotype)]
# levels(cc_mat$y_genotype) <- genotype$short[which(genotype$name %in% cc_mat$y_genotype)]

# # Comparison of proportion of cells being similar and dissimilar
# mod_var <- aov(lm(RA ~ 0 + y_genotype, cc_mat[cc_mat$relation == "variant",]))
# fit_var <- broom::tidy(TukeyHSD(mod_var)) %>%
# arrange(adj.p.value) %>%
# na.omit(.) %>%
# separate(contrast, into = c("x_genotype", "y_genotype"), convert = FALSE, sep = "-") %>%
# # separate(x, into = c("x_genotype", "x_batch"), convert = FALSE, sep = ":") %>%
# # separate(y, into = c("y_genotype", "y_batch"), convert = FALSE, sep = ":") %>%
# # filter(x_batch == y_batch) %>%
# mutate(x_genotype = factor(x_genotype, levels = genotype$short[which(genotype$short %in% .$x_genotype)]),
#   y_genotype = factor(y_genotype, levels = genotype$short[which(genotype$short %in% .$y_genotype)])) %>%
# data.frame(.)

# mod_sim <- aov(lm(RA ~ 0 + y_genotype, cc_mat[cc_mat$relation == "similar",]))
# fit_sim <- broom::tidy(TukeyHSD(mod_sim)) %>%
# arrange(adj.p.value) %>%
# na.omit(.) %>%
# separate(comparison, into = c("x_genotype", "y_genotype"), convert = FALSE, sep = "-") %>%
# # separate(x, into = c("x_genotype", "x_batch"), convert = FALSE, sep = ":") %>%
# # separate(y, into = c("y_genotype", "y_batch"), convert = FALSE, sep = ":") %>%
# # filter(x_batch == y_batch) %>%
# mutate(x_genotype = factor(x_genotype, levels = genotype$short[which(genotype$short %in% .$x_genotype)]),
#   y_genotype = factor(y_genotype, levels = genotype$short[which(genotype$short %in% .$y_genotype)])) %>%
# data.frame(.)

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
                 # facet_wrap(.~ batch, scale = "free", switch = "x") +
                 scale_colour_manual(
                  values = as.character(genotype$`col.idx`), 
                               label = geno.label) +
                 theme_void() +
                  labs(x = "", y = "", colour = "")

  return(ax)


 }, mc.cores = 8)
names(bxp) <- colnames(mds_df)[str_detect(colnames(mds_df), "^mds")]

# Correlational Heatmap
(hmap1 <- cc_samples %>% ggplot(aes(x = source, y = sink)) +
                 geom_raster(aes(fill = saturate(pcc))) +
                 facet_grid(y_genotype ~ x_genotype, space = "fixed", scales = "free", switch = "both") +
                 scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 0.0, 
                                      mid = "white", na.value = "darkgrey") +
                 theme_AKB +
                 theme(axis.text.x = element_blank(),
                       axis.text.y = element_blank(), 
                       strip.text.x = element_text(angle = 90), 
                       strip.text.y = element_text(angle = 180), 
                       legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
                 labs(x = "", y = "", fill = "PCC")) +
                 ggsave(filename = paste(figs, "sample_clustered_heatmap.png", sep = ""), 
                              units = "in",
                           height = 6, 
                           width = 4, device = "png", dpi = 600, 
                        limitsize = FALSE, bg = "transparent")

(mds1_2 <- mds_df %>% 
    ggplot(aes(x= saturate(mds1), 
              y = saturate(mds2))) +
       geom_point(aes(colour = genotype), alpha = 0.6) +
       geom_hline(yintercept = 0.0, alpha = 0.5, 
               colour = "darkgrey", 
               lty = "dashed") +
       geom_vline(xintercept = 0.0, alpha = 0.5, 
               colour = "darkgrey", 
               lty = "dashed") +
       geom_rug(sides = "tr", aes(colour = genotype), position = "identity", alpha = 0.5) +
       # facet_grid(. ~ batch, 
       #  space = "fixed", 
       #  scales = "free", 
       #  switch = "both") +
       scale_colour_manual(values = as.character(genotype$`col.idx`), 
                               label = geno.label) +
       theme_AKB +
       theme() +
       labs(x = paste("MDS-1: ", eigen[1], "%", sep = ""),
            y = paste("MDS-2: ", eigen[2], "%", sep = ""),
            colour = "", 
            size = "")) +
      ggsave(filename = paste(figs, "PQ_sample_MDS12.png", sep = ""), 
                              width = 5, height = 5, 
                              units = "in", device = "png", dpi = 600, 
                        limitsize = FALSE, bg = "transparent")


(mds2_3 <- mds_df %>% 
        ggplot(aes(x= saturate(mds3), 
              y = saturate(mds2))) +
       geom_point(aes(colour = genotype), alpha = 0.6) +
       geom_hline(yintercept = 0.0, alpha = 0.5, 
               colour = "darkgrey", 
               lty = "dashed") +
       geom_vline(xintercept = 0.0, alpha = 0.5, 
               colour = "darkgrey", 
               lty = "dashed") +
       geom_rug(sides = "tr", aes(colour = genotype), position = "identity", alpha = 0.5) +
       # facet_grid(. ~ batch, 
       #  space = "fixed", 
       #  scales = "free", 
       #  switch = "both") +
       scale_colour_manual(values = as.character(genotype$`col.idx`), 
                               label = geno.label) +
       theme_AKB +
       theme() +
       labs(y = paste("MDS-2: ", eigen[2], "%", sep = ""),
            x = paste("MDS-3: ", eigen[3], "%", sep = ""),
            colour = "", 
            size = "")) +
       ggsave(filename = paste(figs, "PQ_sample_MDS23.png", sep = ""), 
                              width = 5, height = 5, 
                              units = "in", device = "png", dpi = 600, 
                        limitsize = FALSE, bg = "transparent")



# Heatmap of features across samples
df_samples_hmap <- cbind.data.frame(SampleID = df_samples$SampleID, 
                                    apply(df_samples[,-1], 2, function(x) (x-mean(x))/sd(x))
                                          ) %>% 
                gather(key = "features", value = "vals", convert = FALSE, -SampleID) %>%
                add_column(genotype = df$genotype[match(.$SampleID, df$SampleID)],
                  batch = df$batch[match(.$SampleID, df$SampleID)]
                  ) %>%
                 mutate(SampleID = factor(SampleID, levels = samples), 
                        features = as.factor(features)
                       ) %>%
                 data.frame(.)

# Sample heatmap, make profile to facet
df_samples_hmap$profile <- ""
df_samples_hmap$profile[str_detect(as.character(df_samples_hmap$feature), "^h")] <- "Haralick"
df_samples_hmap$profile[str_detect(as.character(df_samples_hmap$feature), "^b")] <- "Fluroscence Intensity"
df_samples_hmap$profile[str_detect(as.character(df_samples_hmap$feature), "^s")] <- "Spatial"
                 
(hmap2 <- df_samples_hmap %>% mutate(profile = as.factor(profile)) %>%
ggplot(aes(x = features, y = SampleID)) +
                 geom_raster(aes(fill = saturate(vals))) +
                 facet_grid(genotype ~ profile, space = "free_x", scales = "free", switch = "both") +
                 scale_fill_gradient2(low = "deeppink", high = "darkcyan", na.value = "darkgrey") +
                 theme_AKB +
                 theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5),
                       axis.text.y = element_blank(), 
                       strip.text.x = element_text(angle = 90), 
                       strip.text.y = element_text(angle = 180), 
                       legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
                 labs(x = "", y = "", fill = "Z-score")) +
                 ggsave(filename = paste(figs, "sample_feature_heatmap.png", sep = ""), 
                              width = 4, 
                              height = 6, 
                              device = "png", 
                              units = "in",
                              dpi = 600, 
                        limitsize = FALSE, bg = "transparent")

df_samples <- df_samples %>% 
                gather(key = "features", value = "vals", convert = FALSE, -SampleID) %>%
                cbind.data.frame(., meta[match(.$SampleID, meta$SampleID), which(colnames(meta) %in% facts)]) %>%
                 mutate(SampleID = factor(SampleID, levels = samples), 
                        features = as.factor(features)
                       ) %>%
                 data.frame(.)

levels(df_samples$genotype) <- genotype$short[which(levels(df_samples$genotype) %in% genotype$name)]

# Plot individual samples and their features
stat_obj <- mclapply(unique(as.character(df_samples$features)), function(x){
    
    
#     x <- "h.asm.s1"
    idx <- df_samples$features == x & df_samples$genotype == base
    temp <- df_samples[df_samples$features == x,]
    xmin <- quantile(df_samples$vals[idx], 0.05)
    xmax <- quantile(df_samples$vals[idx], 0.95)
    xmean <- mean(df_samples$vals[idx])
    xmed <- median(df_samples$vals[idx])

    fit <- aov(glm(formula = vals ~ 0 + genotype, data = temp))
    pw.comp <- broom::tidy(TukeyHSD(fit)) %>% mutate(significant = ifelse(adj.p.value < alpha, "*", "NS")) %>% 
    arrange(significant) %>% 
    na.omit(.) %>%
    separate(contrast, into = c("x_genotype", "y_genotype"), convert = FALSE, sep = "-") %>%
    # separate(x, into = c("x_genotype", "x_batch"), convert = FALSE, sep = ":") %>%
    # separate(y, into = c("y_genotype", "y_batch"), convert = FALSE, sep = ":") %>%
    # filter(x_batch == y_batch) %>%
    mutate(term = x) %>%
        data.frame(.)

    return(pw.comp)
    
    
}, mc.cores = 8)

stat_df <- do.call(rbind.data.frame, stat_obj)

# Heatmap for stats
stat_df <- stat_df %>% mutate(significant = as.factor(significant), 
  term = as.factor(term), profile = factor(str_split(term, "\\.", simplify = TRUE)[,1], 
    labels = c("Intensity", "Haralick", "Spatial"), levels = c("b", "h", "s")), 
      logpvals = -log10(adj.p.value), 
      x_genotype = factor(x_genotype, levels = genotype$short[which(genotype$short %in% .$x_genotype)]),
      y_genotype = factor(y_genotype, levels = genotype$short[which(genotype$short %in% .$y_genotype)])
  )

# Build Heatmap for comparative statistics
(hmap3 <- stat_df %>% ggplot(aes(x = x_genotype, y = term, fill = saturate(logpvals))) +
 geom_raster(alpha = 1) +
   geom_tile(fill = NA, aes(colour = significant), width = 0.9, height = 0.9, size = 0.5) +
   facet_grid(profile ~ y_genotype, space = "free", scales = "free", switch = "both") +
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
               )) +
        ggsave(paste(figs, "/Samplewise_features_pvalue_hmap.png", sep = ""), 
               dpi = 600, device = "png", 
               units = "in",
               height = 6, 
               width = 4,  limitsize = F, bg = "transparent")


(hmap4 <- stat_df %>% ggplot(aes(x = x_genotype, y = term, fill = saturate(logpvals))) +
 geom_raster(alpha = 1) +
   geom_tile(fill = NA, aes(colour = significant), width = 0.9, height = 0.9, size = 0.5) +
   facet_grid(profile ~ y_genotype, space = "free", scales = "free", switch = "both") +
   scale_fill_gradient(na.value = "darkgrey", 
     low = "white", 
     high = "red") +
   scale_colour_manual(values = c(`*` = "black", `NS` = "transparent"),
     labels = c(`*` = "FDR < 0.01", `NS` = "")
     ) +
   theme_AKB +
         labs(x = "", 
              y = "", 
              colour = "", fill = "log10(MeanEstimate)") +
         theme(legend.title = element_blank(), 
               axis.title = element_text(hjust = 0.5, vjust = 0.5),
               axis.line = element_blank(),
               axis.text.x = element_text(angle = 90, size = 12, hjust = 0.5, vjust = 0.5),
               strip.text.x = element_text(angle = 90, size = 12, hjust = 0.5, vjust = 0.5),
               axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5)
               )) +
        ggsave(paste(figs, "/Samplewise_features_mean_estimate_hmap.png", sep = ""), 
               dpi = 600, 
               device = "png", 
               units = "in",
               height = 6, 
               width = 4, limitsize = F, bg = "transparent")


# All in one
(all_in_one <- df_samples %>%
    mutate(SampleID = factor(SampleID, levels = samples)) %>%
    ggplot(aes(x = genotype, y = vals)) +
    geom_point(aes(colour = genotype), position = position_jitterdodge(), alpha = 0.1, size = 0.5, shape = 3) + 
    facet_wrap(features ~ ., scales = "free", as.table = TRUE, ncol = 8, nrow = 10) +
    geom_boxplot(fill = NA, 
                 aes(colour = genotype), 
                 alpha = 0.8, 
                 outlier.colour = NA,
                 size = 0.3
                 ) + 
    scale_colour_manual(values = as.character(genotype$col.idx), 
                        label = geno.label) + 
    theme_AKB + 
    theme(axis.text.x = element_blank(),
          strip.text = element_text(size = 8),
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
                         fill = "")) + 
                    ggsave(filename = paste(figs, 
                                            "/sample_feature_summary_plot.png", sep = ""), 
                           width = 8, 
                           height = 7.5, units = "in",
                           device = "png", dpi = 600, 
                           limitsize = FALSE, bg = "transparent")

# Spill canvas
aa <- bxp$mds1 + theme(legend.position = "none") + coord_flip()
bb <- bxp$mds3 + theme(legend.position = "none") + coord_flip()
cc <- bxp$mds2 + theme(legend.position = "none")

# MDS
a1 <- plot_grid(aa, 
  NULL, 
  bb, nrow=1, ncol = 3, 
  rel_height = c(0.1, 0.1, 0.1),
  rel_widths = c(1, 0.2, 1),
  scale = c(0.8)
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
  rel_widths = c(1, 0.4, 1), 
  rel_height = c(.1, 4, .1),
  scale = c(1, 1, 0.5),
  greedy = FALSE
)

ggsave(mds_composite, filename = paste0(figs, "mds_summary.pdf"),
  width = 8, height = 7, 
  units = "in",
  device = "pdf", 
  dpi = 600, 
  limitsize = FALSE, 
  bg = "white")


# Heatmap
a1 <- cowplot::plot_grid(
  hmap1 + theme_void() + theme(legend.position = "none", 
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 180)),
  hmap2 + theme_void() + theme(legend.position = "none", 
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 180)),
  hmap3 + theme_void() + theme(legend.position = "none", 
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 180)),
  hmap4 + theme_void() + theme(legend.position = "none", 
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 180)),
  nrow = 1, ncol = 4, 
  axis = "tblr", 
  rel_widths = c(1, 1, 1, 1), 
  rel_height = c(1, 1, 1)
)

a2 <- cowplot::plot_grid(
  get_legend(hmap1), get_legend(hmap2), get_legend(hmap3), get_legend(hmap3),
  nrow = 1, ncol = 4, 
  align = "hv", 
  rel_widths = c(0.1, .1, .1,.1), 
  rel_height = 1
)

plot_obj <- cowplot::plot_grid(
  a2, a1,
  nrow = 3, ncol = 1, 
  align = "hv", axis = "l",
  scale = c(0.5, 1), 
  greedy = FALSE
)


ggsave(plot_obj, filename = paste0(figs, "profile_summary.pdf"), 
                           width = 20, 
                           height = 20, 
                           device = "pdf", 
                           dpi = 600, 
                           # units = "in",
                           limitsize = FALSE, 
                           bg = "white")

# Each boxplots seperately
if(!dir.exists(paste0(figs, "individual/"))){

  message("Directory created !!")

  dir.create(path = paste0(figs, "individual/"), recursive = TRUE)

} else{

  message("Directory exists !!")
}

mclapply(unique(as.character(df_samples$features)), function(i){

  message(paste0("Spilling ", i))
  
  stat <- c(mn = mean(df_samples$vals[df_samples$features == i]),
            md = median(df_samples$vals[df_samples$features == i]),
            q05 = quantile(df_samples$vals[df_samples$features == i], .05),
            q95 = quantile(df_samples$vals[df_samples$features == i], .95))

  df_samples %>%
      filter(features == i) %>%
      mutate(SampleID = factor(SampleID, levels = samples)) %>%
      ggplot(aes(x = genotype, y = vals)) +
      geom_hline(yintercept = stat, 
        alpha = 0.5, 
                lwd = 0.8,
               colour = "darkgrey", 
               lty = "dashed") +
      geom_point(aes(colour = genotype), position = position_jitterdodge(), alpha = 0.1, size = 0.5, shape = 3) + 
      facet_grid(. ~ batch, scales = "free_x") +
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
                                              "/individual/PQ_sample_(", i,")_feature_plot.png", sep = ""), 
                             width = 6, 
                             height = 4, 
                             units = "in",
                             device = "png", dpi = 600, 
                             limitsize = FALSE, bg = "transparent")

}, mc.cores = 24)
  



message("DONE!!")
sessionInfo()