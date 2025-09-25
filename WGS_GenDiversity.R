# ========== #
#
# Rhipposideros WGS Analysis 2025
#
# Genetic Diversity & Correlations
#
# ========== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(vroom)
library(adegenet)
library(stringr)
library(dplyr)
library(ggplot2)
library(terra)
library(randomcoloR)
library(scales)

# ---------- #
# Genetic Diversity Stats
# ---------- #

# Read in SNP variants
load("data/WGS/snps_hwe_neutral_genlight.RData")

# Calculate observed heterozygosity per individual
Ho <- rowMeans(as.matrix(snps_hwe_neutral) == 1, na.rm = TRUE)
Ho <- Ho[order(names(Ho))] # order individuals alphabetically

# Read in results from ROHan
rohan <- vroom("data/WGS/heterozygosity_ROHan.tsv")

# Subset rohan by the samples that have variants called
diversity_df <- filter(rohan, Ind %in% names(Ho))

# Add site names
diversity_df$Site <- substr(diversity_df$Ind, 1, 3)

# Add variant Ho
diversity_df$Ho_variants <- as.vector(Ho)

# Linear regression between Ho and theta
plot(diversity_df$Het_outside_ROH, diversity_df$Ho_variants)
lm(Het_outside_ROH ~ Ho_variants, data = diversity_df) |> summary()

# Stats per site
(diversity_df_stats <- diversity_df |> 
  group_by(Site) |> 
  summarise(
    Het_outside_ROH_median = median(Het_outside_ROH),
    Het_outside_ROH_mean = mean(Het_outside_ROH),
    Het_inc_ROH_median = median(Het_inc_ROH),
    Het_inc_ROH_mean = mean(Het_inc_ROH),
    Ho_median = median(Ho_variants),
    Ho_mean = mean(Ho_variants),
    Segments_in_ROH_median = median(`Segments_in_ROH_%`),
    Segments_in_ROH_mean = mean(`Segments_in_ROH_%`),
    Avg_length_of_ROH_median = median(Avg_length_of_ROH)
  )
)

# Site order
site_order <- c("PEN","LAN","STO","PAR","BRO","ARL","CAS","HEN","WES","THE","CLA",
                "BAT","ARN","SHE","HID","CAN","PYS","CAD","CRA")

# Change order of sites in data.frame
diversity_df_stats <- mutate(diversity_df_stats, Site = factor(Site, levels = site_order))
diversity_df_stats$Site

# Data for heterozygosity excluding ROHs global medium
global_theta <- data.frame(
  Y = median(diversity_df$Het_outside_ROH),
  Z = "Heterozygosity global median (excluding ROHs)"
)

# Heterozygosity excluding ROHs lollipop graph
(plt_ROH_het <- ggplot(data=diversity_df_stats, aes(x=Site, y=Het_outside_ROH_median))+
  geom_hline(
    data = global_theta,
    aes(yintercept = Y, colour = Z),
    linetype = 2, linewidth = 0.4
  )+
  scale_colour_manual(values = "red")+
  # geom_hline(yintercept=median(rohan$Het_outside_ROH), linetype=3, colour = "red", linewidth=0.5)+
  geom_segment(aes(xend=Site, y=median(diversity_df$Het_outside_ROH), yend=Het_outside_ROH_median), colour="grey")+
  geom_point(size = 3)+
  # geom_point(aes(colour=Site), size=3)+
  # scale_color_manual(values = site_cols)+
  # annotate("text", x = , y = , label = "", size = 3)+
  theme_light()+
  ylab("Heterozygosity (site median excluding ROHs)\n")+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
  )
)

# Save plot
# ggsave(plot = plt_ROH_het, "outputs/WGS_Rhipposideros_heterozygosity_excROHs.png", width = 10, height = 4, dpi = 300)

# Data for heterozygosity excluding ROHs global medium
global_Ho <- data.frame(
  Y = median(diversity_df$Ho_variants),
  Z = "Variant heterozygosity global median"
)

# Variant heterozygosity lollipop graph
(plt_variants_het <- ggplot(data=diversity_df_stats, aes(x=Site, y=Ho_median))+
    geom_hline(
      data = global_Ho,
      aes(yintercept = Y, colour = Z),
      linetype = 2, linewidth = 0.4
    )+
    scale_colour_manual(values = "red")+
    # geom_hline(yintercept=median(rohan$Het_outside_ROH), linetype=3, colour = "red", linewidth=0.5)+
    geom_segment(aes(xend=Site, y=median(diversity_df$Ho_variants), yend=Ho_median), colour="grey")+
    geom_point(size = 3)+
    # geom_point(aes(colour=Site), size=3)+
    # scale_color_manual(values = site_cols)+
    # annotate("text", x = , y = , label = "", size = 3)+
    theme_light()+
    ylab("Variant heterozygosity\n")+
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 10),
    )
)

# Save plot
# ggsave(plot = plt_variants_het, "outputs/WGS_Rhipposideros_heterozygosity_variants.png", width = 10, height = 4, dpi = 300)

# Data for ROH average length
global_ROH_length <- data.frame(
  Y = median(diversity_df$Avg_length_of_ROH),
  Z = "ROH average length global median"
)

# ROHs length lollipop graph
(plt_ROH_length <- ggplot(data=diversity_df_stats, aes(x=Site, y=Avg_length_of_ROH_median))+
    geom_hline(
      data = global_ROH_length,
      aes(yintercept = Y, colour = Z),
      linetype = 2, linewidth = 0.4
    )+
    scale_colour_manual(values = "red")+
    # geom_hline(yintercept=median(rohan$Het_outside_ROH), linetype=3, colour = "red", linewidth=0.5)+
    geom_segment(aes(xend=Site, y=median(diversity_df$Avg_length_of_ROH), yend=Avg_length_of_ROH_median), colour="grey")+
    geom_point(size = 3)+
    # geom_point(aes(colour=Site), size=3)+
    # scale_color_manual(values = site_cols)+
    # annotate("text", x = , y = , label = "", size = 3)+
    theme_light()+
    ylab("ROH average length\n")+
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 10),
    )
)

# Save plot
# ggsave(plot = plt_ROH_length, "outputs/WGS_Rhipposideros_ROH_length.png", width = 10, height = 4, dpi = 300)

# ROHs % lollipop graph
(plt_ROH_length <- ggplot(data=diversity_df_stats, aes(x=Site, y=Segments_in_ROH_median))+
    # geom_hline(
    #   data = global_ROH_length,
    #   aes(yintercept = Y, colour = Z),
    #   linetype = 2, linewidth = 0.4
    # )+
    scale_colour_manual(values = "red")+
    # geom_hline(yintercept=median(rohan$Het_outside_ROH), linetype=3, colour = "red", linewidth=0.5)+
    geom_segment(aes(xend=Site, y=median(diversity_df$`Segments_in_ROH_%`), yend=Segments_in_ROH_median), colour="grey")+
    geom_point(size = 3)+
    # geom_point(aes(colour=Site), size=3)+
    # scale_color_manual(values = site_cols)+
    # annotate("text", x = , y = , label = "", size = 3)+
    theme_light()+
    ylab("ROH average length\n")+
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 10),
    )
)

# ---------- #
# Prepare Environment Data
# ---------- #

# Read in coordinates in British National Grid projection
coords <- read.csv("../Coordinates/coordinates_BritishNationalGrid.csv")
coords_vect <- vect(coords, geom=c("Lon","Lat"), crs="epsg:27700")
plot(coords_vect)

# Light raster
nightlight <- rast("data/WGS/spatial_data/light_tiff.tif")
nightlight <- project(nightlight, "epsg:27700")
# nightlight; plot(nightlight)

# Create a 2km buffer
coords_buffer <- buffer(coords_vect, width=2000)

# Create data.frame to store results
envir_df <- data.frame(Site = coords$Site_old)

# Extract the mean light within 2km buffer for each site
nightlight_vals <- terra::extract(nightlight, coords_buffer, fun = mean, na.rm = TRUE)
envir_df$nightlight <- nightlight_vals$light_tiff

# Read in 1km Land Cover 2023 GeoTIFF
land_cover <- rast("data/WGS/spatial_data/gb2023lcm1km_percentage_aggregate.tif")

# Classes to extract
classes <- list(
  "broadleaf_woodland" = 1,
  # "coniferous_woodland" = 2,
  "all_woodland" = c(1,2),
  "arable" = 3,
  "urban" = 10
)

# Function to extract values within 2km buffer from land cover classes
extract_class_values <- function (raster, class, coordinates) {
  values_df <- terra::extract(raster[[class]], coordinates, fun = mean, na.rm = TRUE)
  
  if (length(class) == 1) {
    return(values_df[[2]])  
  } else {
    return(
      values_df |> 
        dplyr::select(starts_with("gb2023lcm1km_percentage_aggregate")) |> 
        mutate(sums = rowSums(across(everything()))) |> 
        pull(sums)
    )
  }
}

# Extract values for classes
land_cover_vals <- classes |> 
  lapply(X = _, \(x) extract_class_values(land_cover, x, coords_buffer)) |> 
  list2DF()

# Add land cover data to data.frame
envir_df_lc <- cbind(envir_df, land_cover_vals)
head(envir_df_lc)

# Read in CHELSA climate rasters
chelsa_filepaths <- list.files("data/WGS/spatial_data/", pattern = "UK_CHELSA", full.names = T)
chelsa_rasters <- lapply(chelsa_filepaths, rast)
names(chelsa_rasters) <- str_extract(chelsa_filepaths, "(?<=UK_CHELSA_)[^_]+")

# Extract values for climate rasters
chelsa_vals <- sapply(names(chelsa_rasters), \(x) terra::extract(chelsa_rasters[[x]], project(coords_vect, "EPSG:4326"))[[2]] )

# Read in National Forest Inventory 2023 vector GeoPackage
forest <- vect("data/WGS/spatial_data/National_Forest_Inventory_GB_2023_study_coords_5km_buffer.gpkg")
# plot(forest)

# Find the nearest woodland feature to each coordinate
forest_dist <- nearest(coords_vect, forest)
forest_dist$distance

# Add CHELSA data to data.frame
envir_df_final <- cbind(envir_df_lc, woodland_dist = forest_dist$distance, chelsa_vals)
head(envir_df_final)

# Transform environmental and climate data using log1p
# envir_df_log <- as.data.frame(apply(envir_df[,2:ncol(envir_df)], 2, log1p))
# colnames(envir_df_log) <- paste0(colnames(envir_df_log), "_log")
# envir_df_log$Site <- envir_df$Site
# envir_df_log

# ---------- #
# Prepare Data Frames
# ---------- #

# Prepare site-based data.frame for correlation tests
df <- data.frame(
  envir_df_final,
  # envir_df[-1],
  # envir_df_log,
  Heterozygosity = diversity_df_stats$Het_outside_ROH_median,
  Ho_variants = diversity_df_stats$Ho_median,
  ROH_length =  diversity_df_stats$Avg_length_of_ROH_median,
  ROH_segments = diversity_df_stats$Segments_in_ROH_median
)
df

# Read in coordinates, add to data.frame and export CSV
read.csv("../Coordinates/coordinates.csv") |> 
  left_join(df, by = join_by("Site_old" == "Site")) |> 
  write.csv(file = "outputs/WGS_Rhipposideros_heterozygosity_site.csv", row.names = F)

# Prepare individual-based data.frame for correlation tests
df <- left_join(diversity_df, envir_df_final, by = "Site") |> 
  # left_join(envir_df_final, by = "Site") |> 
  # Remove BRO10 outlier
  filter(Ind != "BRO10")
df

# Read in coordinates, add to data.frame and export CSV
read.csv("../Coordinates/coordinates.csv") |> 
  left_join(df, by = join_by("Site_old" == "Site")) |> 
  write.csv(file = "outputs/WGS_Rhipposideros_heterozygosity_ind.csv", row.names = F)

# ---------- #
# Linear Regression
# ---------- #

# Theta heterozygosity with nightlight
plot(df$Het_outside_ROH, df$nightlight)
cor.test(df$Het_outside_ROH, df$nightlight, method = "spearman")

# Variant heterozygosity with nightlight
plot(df$Ho_variants, df$nightlight)
cor.test(df$Ho_variants, df$nightlight, method = "spearman")

# ROH length with nightlight
plot(df$Avg_length_of_ROH, df$nightlight)
cor.test(df$Avg_length_of_ROH, df$nightlight, method = "spearman")

# ROH % with nightlight
plot(df$`Segments_in_ROH_%`, df$nightlight)
cor.test(df$`Segments_in_ROH_%`, df$nightlight, method = "spearman")

# Theta heterozygosity with broadleaf woodland
plot(df$Het_outside_ROH, df$broadleaf_woodland)
cor.test(df$Het_outside_ROH, df$broadleaf_woodland, method = "spearman")

# Variant heterozygosity with broadleaf woodland
plot(df$Ho_variants, df$broadleaf_woodland)
cor.test(df$Ho_variants, df$broadleaf_woodland, method = "spearman")

# ROH length with broadleaf woodland
plot(df$Avg_length_of_ROH, df$broadleaf_woodland)
cor.test(df$Avg_length_of_ROH, df$broadleaf_woodland, method = "spearman")

# ROH % with broadleaf woodland
plot(df$`Segments_in_ROH_%`, df$broadleaf_woodland)
cor.test(df$`Segments_in_ROH_%`, df$broadleaf_woodland, method = "spearman")

# Theta heterozygosity with all woodland
plot(df$Het_outside_ROH, df$all_woodland)
cor.test(df$Het_outside_ROH, df$all_woodland, method = "spearman")

# Variant heterozygosity with all woodland
plot(df$Ho_variants, df$all_woodland)
cor.test(df$Ho_variants, df$all_woodland, method = "spearman")

# ROH length with all woodland
plot(df$Avg_length_of_ROH, df$all_woodland)
cor.test(df$Avg_length_of_ROH, df$all_woodland, method = "spearman")

# ROH % with all woodland
plot(df$`Segments_in_ROH_%`, df$all_woodland)
cor.test(df$`Segments_in_ROH_%`, df$all_woodland, method = "spearman")

# Theta heterozygosity with distance to woodland
plot(df$Het_outside_ROH, df$woodland_dist)
cor.test(df$Het_outside_ROH, df$woodland_dist, method = "spearman")

# Variant heterozygosity with distance to woodland
plot(df$Ho_variants, df$woodland_dist)
cor.test(df$Ho_variants, df$woodland_dist, method = "spearman")

# ROH length with distance to woodland
plot(df$Avg_length_of_ROH, df$woodland_dist)
cor.test(df$Avg_length_of_ROH, df$woodland_dist, method = "spearman")

# ROH % with distance to woodland
plot(df$`Segments_in_ROH_%`, df$woodland_dist)
cor.test(df$`Segments_in_ROH_%`, df$woodland_dist, method = "spearman")

# Theta heterozygosity with arable
plot(df$Het_outside_ROH, df$arable)
cor.test(df$Het_outside_ROH, df$arable, method = "spearman")

# Variant heterozygosity with arable
plot(df$Ho_variants, df$arable)
cor.test(df$Ho_variants, df$arable, method = "spearman")

# ROH length with arable
plot(df$Avg_length_of_ROH, df$arable)
cor.test(df$Avg_length_of_ROH, df$arable, method = "spearman")

# ROH % with arable
plot(df$`Segments_in_ROH_%`, df$arable)
cor.test(df$`Segments_in_ROH_%`, df$arable, method = "spearman")

# Theta heterozygosity with urban
plot(df$Het_outside_ROH, df$urban)
cor.test(df$Het_outside_ROH, df$urban, method = "spearman")

# Variant heterozygosity with urban
plot(df$Ho_variants, df$urban)
cor.test(df$Ho_variants, df$urban, method = "spearman")

# ROH length with urban
plot(df$Avg_length_of_ROH, df$urban)
cor.test(df$Avg_length_of_ROH, df$urban, method = "spearman")

# ROH % with urban
plot(df$`Segments_in_ROH_%`, df$urban)
cor.test(df$`Segments_in_ROH_%`, df$urban, method = "spearman")
