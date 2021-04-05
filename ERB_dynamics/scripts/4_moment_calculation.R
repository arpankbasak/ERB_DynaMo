rm(list = ls())

# Run this script for annotating the grayscale images only
# For colour donot use thresholding, instead treat the red channel as potential cellwall eroding the lines
message(paste("---", "==Computing Feature Statistics=="," Job Started at ", as.character(Sys.time()), "---", sep = ""))

# Load Requisites
pkgs <- c("EBImage", "tidyverse", "parallel", "vegan", "cowplot", "gridExtra", "patchwork")
lapply(pkgs, require, character.only = T)
options(warn = 1, mc.cores = 8)

# Measures not to fry R
options(set.cores = 8)
param <- "./script/parameters.R"
path <- as.character(getwd())
setwd(path)
source(param)

figs <- paste0(figs, "sample_wise_analysis_")

load(paste(data, "Obj_image_features.RData", sep = ""))

mmat <- as.matrix(imat[,grep(colnames(imat), pattern = "^m")])

# Read metadata
meta <- read.table(paste(output, "metadata.txt", sep = ""), sep = "\t", header = T, as.is = T)

cname <- facts

# meta$genotype <- str_split(meta$SampleID, pattern = "_", simplify = T)[,2]
id <- which(!str_detect(row.names(mmat), "\\."))
row.names(mmat)[id] <- paste(row.names(mmat)[id], "0", sep = ".")

# Make the dataframe
df <- cbind.data.frame(FeatureID = row.names(mmat), mmat) %>%
    mutate(FeatureID = str_replace_all(FeatureID, "feature-", "")) %>%
    separate(FeatureID, into = c("SampleID", "feature_no"), sep = "\\.", remove = FALSE, convert = FALSE) %>%
    separate(SampleID, remove = FALSE, into = facts, sep = "_") %>%
    mutate(m.cx = as.numeric(m.cx), 
      m.cy = as.numeric(m.cy),
      m.majoraxis = as.numeric(m.majoraxis),
      m.eccentricity = as.numeric(m.eccentricity),
      m.theta = as.numeric(m.theta)) %>%
    data.frame(.)

df <- df[!str_detect(df$SampleID, "_root_"),]

# Extract out morphology features eccentricity
levels(meta$genotype) <- as.character(genotype$name)
idx <- which(genotype$short %in% df$genotype)
geno.label <- geno.label[idx]
genotype <- genotype[idx,]
geno <- genotype

df$genotype <- factor(df$genotype, levels = genotype$short)
df$batch <- as.factor(df$batch)
meta$genotype <- factor(meta$genotype, levels = genotype$short)
df$SampleID <- as.factor(df$SampleID)
df$batch <- as.factor(df$batch)
df$depth <- as.factor(df$depth)
df$time <- as.factor(df$time)


# Momentary calculation
df_m <- df %>% mutate(cm = abs(sqrt((m.cx^2) + (m.cy^2))/(m.cx + m.cy))) %>%
separate(SampleID, into = c("SampleID_t", "x"), remove = FALSE, sep = "_t") %>%
mutate(SampleID_tn = paste(SampleID_t, feature_no, sep = "_")) %>%
select(-x)

df_disp <- df_m %>% 
select(SampleID_tn, time, cm) %>%
spread(key = time, value = cm, fill = 0) %>%
data.frame(.)

df_m$genotype <- as.character(meta$genotype)[match(as.character(df_m$SampleID), as.character(meta$SampleID))]

df_m$genotype <- factor(df_m$genotype, levels = genotype$short[genotype$short %in% df_m$genotype])

row.names(df_disp) <- df_disp[,1]
vel <- rowSums(df_disp[,-1])/10
temp <- apply(df_disp[,-c(1)], 1, function(x) (x-x[1]))
temp <- apply(t(temp), 2, function(x) ma(x, length(x))) # moving average
temp[,"t0"] <- 0

time_df <- cbind.data.frame(id = row.names(temp), temp) %>%
gather(key = "time", value = "vals", -id) %>%
separate(id, remove = FALSE, into = c(facts[-length(facts)], "no_feature"), sep = "_") %>%
mutate(genotype = factor(.$genotype, levels = as.character(geno$short)),
  replicate = as.factor(replicate),
  depth = as.factor(depth),
  batch = as.factor(batch),
  time = as.numeric(str_replace_all(time, "t", "")),
  velocity = abs(vals)/10,
  vals = abs(vals)) %>%
data.frame(.)

mean_time_df <- time_df %>% group_by(genotype, depth, time) %>%
summarise(
  mean_disp = (mean(vals))^2, median_disp = (median(vals))^2, 
  mean_vel = (mean(velocity)), median_vel = (median(velocity)), 
  acc = (mean(velocity))^2) %>%
data.frame(.)

# boxplot
time_df %>%
mutate(zer = factor(as.character(time), levels = rev(0:9))) %>%
ggplot(aes(x = genotype, y = abs(vals))) +
geom_hline(yintercept = 0, colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_point(size = 0.5, aes(colour = genotype), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.4) +
geom_boxplot(aes(colour = genotype), fill = NA, width = 0.8, outlier.alpha = 0) +
facet_grid(zer~., 
        space = "free_x", 
        scales = "free", 
        switch = "both") +
scale_colour_manual(values = as.character(genotype$col.idx), 
                       label = geno.label) +
theme_AKB +
theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5),
  axis.text.y = element_text(size = 10)) +
labs(x = "",
    y = "Displacement",
    colour = "", 
    size = "") +
ggsave(filename = paste("./figures/displacement_boxplot.png", sep = ""), 
                      width = 3, height = 7, 
      device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")

# Line plot -- Displacement
time_df %>%
mutate(random = as.factor(paste(replicate, sep = "_"))) %>%
ggplot(aes(x = time, y = vals^2, colour = genotype)) +
geom_hline(yintercept = 0, colour = "darkgrey", lty = "solid", lwd = 0.8) +
# geom_line(aes(colour = genotype, group = random), lty = "solid", alpha = 0.4) +
geom_point(size = 0.5, aes(fill = genotype), colour = "black", alpha = 0.6) +
geom_smooth(data = mean_time_df, aes(x = time, y = median_disp, colour = genotype), method = "loess", se = FALSE) +
geom_rug(sides = "r", aes(colour = genotype), position = "identity", alpha = 0.5) +
# facet_grid(depth ~ ., 
#         space = "free_x", 
#         scales = "free", 
#         switch = "both") +
scale_colour_manual(values = as.character(genotype$col.idx), 
                       label = geno.label) +
scale_fill_manual(values = as.character(genotype$col.idx), 
                       label = geno.label) +
theme_AKB +
theme() +
labs(x = "time",
    y = "Moving average displacement",
    colour = "", 
    size = "") +
ggsave(filename = paste("./figures/moving_av.png", sep = ""), 
                      width = 5, height = 3, 
      device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")

# Velovity
time_df %>%
# mutate(random = as.factor(paste(depth, replicate, sep = "_"))) %>%
ggplot(aes(x = time, y = velocity)) +
geom_hline(yintercept = 0, colour = "darkgrey", lty = "solid", lwd = 0.8) +
# geom_line(aes(colour = genotype, group = random), lty = "solid", alpha = 0.4) +
geom_point(size = 0.5, aes(colour = genotype), alpha = 0.6) +
geom_smooth(data = mean_time_df, aes(x = time, y = mean_vel, colour = genotype), method = "loess", se = FALSE) +
geom_rug(sides = "r", aes(colour = genotype), position = "identity", alpha = 0.5) +
# facet_grid(depth ~ ., 
#         space = "free_x", 
#         scales = "free", 
#         switch = "both") +
scale_colour_manual(values = as.character(genotype$col.idx), 
                       label = geno.label) +
theme_AKB +
theme() +
labs(x = "time",
    y = "Velocity",
    colour = "", 
    size = "") +
ggsave(filename = paste("./figures/moving_av_velocity.png", sep = ""), 
                      width = 5, height = 3, 
      device = "png", dpi = 600, limitsize = FALSE, bg = "transparent")


# Heatmap for displacement
temp_df <- cbind.data.frame(id = row.names(temp), temp) %>%
gather(key = "time", value = "vals", -id) %>%
separate(id, remove = FALSE, into = c(facts[-length(facts)], "no_feature"), sep = "_") %>%
mutate(genotype = as.factor(genotype),
  replicate = as.factor(replicate),
  depth = as.factor(depth),
  batch = as.factor(batch),
  time = as.factor(time)) %>%
data.frame(.)

# Plot the stats -- try loess regression or something of time series
temp_gen_df <- time_df %>% group_by(genotype, replicate, depth, time) %>%
summarise(vals = sum(abs(vals))) %>%
ungroup() %>%
mutate(time = as.factor(as.character(time))) %>%
data.frame(.)

mod_glm_cum <- glm(log2(vals+1) ~ 0 + genotype+time, data = temp_gen_df)
# mod_loess_cum <- loess(log2(vals+1) ~ 0 + time, data = temp_gen_df, parametric = TRUE)
mod <- aov(mod_glm_cum)
fit <- broom::tidy(TukeyHSD(mod)) %>% 
separate(contrast, into = c("A", "B")) %>% 
mutate(logpvals = -log10(adj.p.value), 
  significance = adj.p.value < 0.05,
  A = as.factor(A),
  B = as.factor(B)) %>%
arrange(A,B, significance) %>%
data.frame(.)

geno <- genotype
temp_gen_df %>%
mutate(time = as.numeric(str_replace_all(as.character(time), "t", "")),
  genotype = factor(genotype, levels = geno$short)) %>%
ggplot(aes(x= time, y = abs(vals), fill = genotype, colour = genotype)) +
geom_point(size = 1, shape = 23, colour = "black") +
geom_smooth(method = "loess", se = FALSE) +
# facet_grid(depth ~ ., space = "free", scales = "fixed", switch = "both") +
scale_fill_manual(values = geno$col.idx) +
scale_colour_manual(values = geno$col.idx) +
theme_AKB +
theme() +
labs(x = "time", 
  y = "displacement", 
  colour = "") +
                 ggsave(filename = paste("./figures/trendline_loess.png", sep = ""), 
                              width = 7, height = 4, device = "png", dpi = 600, unit = "in",
                        limitsize = FALSE, bg = "transparent")

(hmap_p <- fit %>%
  filter(term == "genotype") %>%
  mutate(A = factor(A, levels = genotype$short),
    B = factor(B, levels = genotype$short)) %>% 
  ggplot(aes(x = A, y = B)) +
  geom_raster(aes(fill = (logpvals))) +
  geom_tile(aes(color = significance), width = 0.9, height = 0.9, size = 1, fill = NA) +
 # facet_grid(genotype ~ ., space = "free_x", scales = "free", switch = "both") +
 scale_fill_gradient2(low = "lightgrey", midpoint = -log10(.01), 
  mid = "white", 
  high = "darkred", na.value = "white") +
 scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA)) +
 theme_AKB +
 theme(axis.text.x = element_text(size = 12),
       axis.text.y = element_text(size = 12), 
       strip.text.x = element_text(angle = 90), 
       strip.text.y = element_text(angle = 180), 
       legend.text = element_text(hjust = 0.5, vjust = 0.5)) +
 labs(x = "", y = "", fill = "")) +
 ggsave(filename = paste("./figures/stats_glm.png", sep = ""), 
              width = 4, 
              height = 4, 
              units = "in",
              device = "png", dpi = 600, 
        limitsize = FALSE, bg = "transparent")



message("DONE!!")
sessionInfo()