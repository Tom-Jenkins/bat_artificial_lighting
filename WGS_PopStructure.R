# ========== #
#
# Rhipposideros WGS Analysis 2025
#
# Population Structure
#
# ========== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(adegenet)
library(dartR)
# dartR::gl.install.vanilla.dartR() # install dartR dependencies
library(LEA)
library(mapmixture)
library(ggspatial)
library(sf)
library(terra)
library(stringr)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthhires)
library(gridExtra)

# ---------- #
# Principal Components Analysis
# ---------- #

# Read in SNPs
load("data/WGS/snps_hwe_neutral_genlight.RData")
snps_hwe_neutral

# PCA
x = tab(snps_hwe_neutral, NA.method = "mean")
pca1 = dudi.pca(x, scannf = FALSE, scale = TRUE, nf = 3)

# Extract percent of genetic variance is explained by each axis
(percent <- round(pca1$eig/sum(pca1$eig)*100, 1))

# Create a divergent colour palette
# set.seed(1234)
# library(randomcoloR)
# cols <- distinctColorPalette(nPop(snps_hwe_neutral))
cols <- c(
  "#E3A6A4", "#E05389", "#C1EB4E", "#797BDC", "#D3D6DE", "#6EE165",
  "#7F7A8F", "#97936A", "#D654D3", "#DEC15A", "#D482D2", "#70E2AE",
  "#963CE6", "#E07C4B", "#E0E1B5", "#D2B1E6", "#72B2DC", "#C0DF88",
  "#86DCD8"
)

# Change site labels to new labelling scheme for paper
codes <- c(
  "DE1", "BR", "SO1", "DE2", "NW1", "HE", "SO2", "SO3", 
  "NW2", "SO4", "GL1", "CO1", "DE3", "CO2", "SW", "GL2", 
  "CO3", "BE", "SO5"
)
levels(snps_hwe_neutral$pop) <- codes
site_order <- c("CO2","CO1","CO3","DE3","DE2","DE1","SO2","SO4","SO5","BE","SO3",
                "SO1","BR","GL1","GL2","HE","SW","NW1","NW2")
snps_hwe_neutral$pop <- factor(snps_hwe_neutral$pop, levels = site_order)
snps_hwe_neutral$pop

# Data frame
df <- data.frame(Sample = indNames(snps_hwe_neutral), Site = snps_hwe_neutral$pop, pca1$li)
head(df)

# Plot
# scatter_plot(
#   dataframe = df[,3:ncol(df)],
#   group_ids = snps_hwe_neutral$pop,
#   colours = cols, point_size = 5,
#   type = "labels", labels = snps_hwe_neutral$pop, segments = F, centroids = F,
#   plot_title = paste("Lesser horseshoe bats PCA: ", snps_hwe_neutral$n.loc,"biallelic SNPs"),
#   percent = percent
# )
(plt_pca <- scatter_plot(
  dataframe = df[,3:ncol(df)],
  group_ids = snps_hwe_neutral$pop,
  colours = cols, point_size = 5, centroid_size = 3,
  type = "points", labels = codes, segments = F, centroids = T,
  # plot_title = paste("Lesser horseshoe bats PCA: ", snps_hwe_neutral$n.loc,"biallelic SNPs"),
  percent = percent, axes = c(1,2)
)+
  theme(
    plot.title = element_text(size = 15),
  )+
  annotate("text", x = -45, y = -70, label = "North Wales", size = 6)+
  # annotate("text", x = -60, y = 10, label = "South Wales", size = 6)+
  annotate("text", x = 30, y = 30, label = "England", size = 6)
)

# Save plot
ggsave(plot = plt_pca, "outputs/WGS_Rhipposideros_PCA.png", width = 10, height = 6, dpi = 300)

# ---------- #
# SNMF
# ---------- #

# If true re-run snmf
runSNMF <- FALSE

# Run snmf algorithm
if (runSNMF) {
  snmf1 <- snmf(
    input.file = "data/WGS/snps_hwe_neutral.geno",
    K = 1:10, # number of K ancestral populations to run
    repetitions = 10, # ten repetitions for each K
    entropy = TRUE, # calculate cross-entropy
    project = "new",
    ploidy = 2,
    seed = 1234,
    CPU = 4
  )
}

# Load snmf project
snmf1 <- load.snmfProject("data/WGS/snps_hwe_neutral.snmfProject")

# Plot cross-entropy results to assess optimal number of K
# Smaller values of cross-entropy usually mean better runs
# A plateau usually represents the K that best fits the data
plot(snmf1, col = "blue", cex = 1.5, pch = 19, cex.axis = 1, cex.lab = 1.5,
     main = "SNMF: Cross-entropy", cex.main = 1.5)

# If true export entropy plot
exportEntropy <- FALSE
if (exportEntropy) {
  png("outputs/WGS_Rhipposideros_entropy.png", width = 10, height = 7, unit = "in", res = 300)
  plot(snmf1, col = "blue", cex = 1.5, pch = 19, cex.axis = 1, cex.lab = 1.5,
       main = "SNMF: Cross-entropy", cex.main = 1.5)
  dev.off()
  pdf("outputs/WGS_Rhipposideros_entropy.pdf", width = 10, height = 7)
  plot(snmf1, col = "blue", cex = 1.5, pch = 19, cex.axis = 1, cex.lab = 1.5,
       main = "SNMF: Cross-entropy", cex.main = 1.5)
  dev.off()
}

# Extract Q-matrix for each K using best run
best_K <- sapply(1:10, \(i) which.min(cross.entropy(snmf1, K = i)))
qmatrix_ls <- lapply(1:10, \(x) Q(snmf1, K = x, run = best_K[x]))
names(qmatrix_ls) <- paste0("K", 1:10)

# Create data.frame to store admixture results
snmf_df <- data.frame(
  Site = snps_hwe_neutral$pop,
  Ind = snps_hwe_neutral$ind.names
)

# K to plot
K <- 3

# Add admixture results to new data.frame
snmf_admix <- cbind(snmf_df, qmatrix_ls[[paste0("K",K)]]) |> arrange(Site)
head(snmf_admix)

# Cluster colours
k_cols <- c("#fc8d62","#66c2a5","#8da0cb")

# Plot
(plt_structure <- structure_plot(
  admixture_df = snmf_admix,
  cluster_cols = k_cols,
  site_order = site_order,
  site_labels_size = 4,
  site_ticks = FALSE,
  site_dividers = TRUE, divider_width = 0.5,
  site_labels_y = -0.12,
  ylabel = "Ancestry proportion"
))

# Save plot
ggsave(plot = plt_structure, "outputs/WGS_Rhipposideros_structure.png", width = 10, height = 4, dpi = 300)

# Read in coordinates and add columns for transformed coordinates
coords <- read.csv("../Coordinates/coordinates.csv")
head(coords)

# Add columns for transformed coordinates
coords <- st_as_sf(coords, coords = c("Lon","Lat"), crs = 4326) |> 
  st_transform(crs = 27700) |> 
  st_coordinates() |> 
  cbind(coords)
head(coords)

# Adjust coordinates for plotting
coords <- coords %>%
  mutate(
    X_adjust = case_when(
      Site == "BR" ~ X + 18000,
      Site == "SO1" ~ X + 19000,
      Site == "BE" ~ X + 20000,
      .default = X
    ),
    Y_adjust = case_when(
      Site == "BR" ~ Y + 13000,
      Site == "SO1" ~ Y + 8000,
      Site == "BE" ~ Y,
      .default = Y + 14000
    )
  )
head(coords)

# High resolution countries boundary vector layer
basemap <- rnaturalearthhires::countries10[,"geometry"]

# Mapmixture
(plt_map <- mapmixture(
  admixture_df = snmf_admix,
  coords_df = coords[, c("Site","Lat","Lon")],
  basemap = basemap,
  crs = 27700,
  boundary = c(xmin=-5.75, xmax=-1.25, ymin=49.9, ymax=53.5),
  pie_size = 0.15,
  cluster_cols = k_cols,
  cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
  # sea_colour = "#deebf7",
  # land_colour = "#d9d9d9",
)+
  theme(
    axis.title = element_text(size = 12),
  )+
  geom_label(
    data = coords,
    aes(x = X_adjust, y = Y_adjust, label = Site),
    size = 3
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 6, alpha = 1)))
)

# Save plot
ggsave(plot = plt_map, "outputs/WGS_Rhipposideros_mapmixture.png", width = 8, height = 10, dpi = 300)

# ---------- #
# Figure
# ---------- #

# Layout matrix
layout_mat <- rbind(
  c(1,2),
  c(3,3)
)

# Arrange plots
plt_figure <- grid.arrange(
  plt_pca+ labs(tag = "A"),
  plt_map+ labs(tag = "B"),
  plt_structure+ labs(tag = "C"),
  ncol = 2, nrow = 2, layout_matrix = layout_mat, heights = c(2,1)
)

# Save figure
ggsave(
  plot = plt_figure,
  filename = "figures/Figure3.png",
  width = 15, height = 10, dpi = 600
)
ggsave(
  plot = plt_figure,
  filename = "figures/Figure3.pdf",
  width = 15, height = 10
)
