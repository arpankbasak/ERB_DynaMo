# Set path
data = "./data/"
image = "./images/"
output = "./output/"
# path = "/klaster/work/abasak/git_repo_erb/ERB_ImageAnalysis/ERB_dynamics/"
stats = "./statistics/"
figs = "./figures/"
cell_image = "./images/segmented_images/stacked/cells/"
feature_image = "./images/segmented_images/stacked/features/"
live_stock = "./images/live_stock/"
live_stock = "./images/merged_images/"
feature_image = "./images/clustered_images/"
dyn_image = "./images/dynamics/"
cell_dyn_image = "./images/dynamics/cells/"
feature_mat = "./output/feature_matrix/"


# Prerequisites
facts <- c("date", "genotype", "tissue_type", "dop", "batch", "replicate", "depth", "time")
feature_type <- c("erb", "cell")
feature_names <- c("haralick_features", "spatial_features", "moment_features", "intensity_features")
groups_genotype <- c("nai1-gfph", "gfph")
eval_genotype <- c("bglu21-1-leb-1GFPHF22", "meb2-1", "meb1-2", "meb1-1meb2-1")

col_hex <- c("#8B008B",  "#006400", "#BD3430", "#4A2D4E","#875525","#A83683", "#4E655E", "#853541","#3A3120", "#535D8E", "#476B2A", "#7851A8", "#BD3230", "#4A2A4E","#875225","#A83783", "#4E635E", "#853241","#3A3320", "#535D9E", "#476C2A", "#7851B8", "#BD3330", "#4A2B4E","#875325","#A83883", "#4E625E", "#853341","#3A3420", "#535D7E", "#476D2A", "#7851C8", "#BD3130", "#4A2C4E","#875425","#A83983", "#4E615E", "#853141","#3A3220", "#535D6E")

# Table for the statistical parameters
seed <- 1
k_folds <- 25
epoch <- 100
batch_size <- 100
val_split <- 0.1
n_factors <- 32
px_scale <- 0.42
noise <- 1000
p.adjust.method <- "fdr"
partition <- 0.05
alpha <- 0.01
base <- "nai1"

# Table for genotype
genotype <- data.frame(name = c("nai1-gfph", "nai1-pPYK10-tdTOM-GFPh", "gfph", "GFPh", 
                                "bglu21-1-leb-1GFPHF22", "meb2-1", "meb1-2", "meb1-1meb2-1",
                       "lgo-2GFPhF31", "dek1"),
                      short = c("nai1", "nai1.pPYK10", "gfph", "GFPh", 
                                "leb1", "meb2", "meb1", "meb1meb2",
                       "lgo2", "dek1"),
                      shorter = c("nai1", "nai1", "gfph", "gfph", 
                                "leb1", "meb2", "meb1", "meb1meb2",
                       "lgo2", "dek1"),
                      com = c("No ERBs", "No ERBs", "Wild type", "Wild type", 
                                "Long ERBs", "Irregular ERBs", "Irregular ERBs", "Irregular ERBs",
                       "Defect in cell", "Defect in cell"),
                       col.idx = c("darkgreen", "palegreen", "darkmagenta", "magenta",
                                   "darkcyan", "darkkhaki", "khaki", "greenyellow", "darkorange", "turquoise"), 
                       shape.idx = c(1, 2, 3, 4, 
                                     6, 8, 11, 12, 13, 14), stringsAsFactors = FALSE)

# Table for batch
batch <- data.frame(name = c("Experiment1", "Experiment2"),
  col.idx = c("lightgrey", "black"),
  shape.idx = c(2, 6),
  description = c("","")
  )

# Fancy genotype
geno.label <- c(expression(italic("nai1")),
                expression(italic("nai1-pPYK10")),
                expression(italic("gfph")),
                expression(italic("GFPh")),
                expression(italic("bglu21-1-leb-1GFPHF22")),
                expression(italic("meb2-1")),
                expression(italic("meb1-2")),
                expression(italic("meb1-1meb2-1")),
                expression(italic("lgo2")),
                expression(italic("dek1"))
)

# Featur
features <- data.frame(name = c("cell", "erb"), 
                       col.idx = c("lightred", "lightblue"), 
                       shape.idx = c(3, 5), stringAsFactors = FALSE)

# Theme set
theme_AKB <- theme_bw() +
    theme(text = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          legend.key=element_blank(),
          legend.text.align=0,
          legend.text=element_text(size = 8, hjust = 0.5, vjust = 0.5),
          legend.position="top",
          strip.text=element_text(face="bold", size = 16),
          axis.line=element_line(),
          axis.line.x=element_line(size = 1),
          axis.line.y=element_line(size = 1),
          panel.background=element_rect(fill="transparent", colour=NA),
          plot.background=element_rect(fill="transparent", colour=NA),
          strip.background=element_rect(fill="transparent", colour=NA),
          strip.placement="outside")

# Functions
# Saturation of intensities especially for plotting heatmaps
saturate <- function(x){
    max <- quantile(x, .99)
    min <- quantile(x, .01)
    
    idx <- (x > max)
    x[idx] <- max
    
    idx <- x < min
    x[idx] <- min
    return(x)
}

# Optimum Clusters
kmeansAIC <- function(fit){
    m <- ncol(fit$centers)
    n <- length(fit$cluster)
    k <- nrow(fit$centers)
    D <- fit$tot.withinss
    return(c(
        AIC=(D + 2*m*k),
        BIC=(D + log(n)*m*k)
        )
    )
}

ma <- function(arr, n=15){
  res = arr
  for(i in n:length(arr)){
    res[i] = mean(arr[(i-n+1):i])
  }
  return(res)
}

