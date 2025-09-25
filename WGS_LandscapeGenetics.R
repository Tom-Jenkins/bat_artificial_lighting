# ========== #
#
# Rhipposideros WGS Analysis 2025
#
# Landscape Genetics
#
# ========== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(adegenet)
library(ResistanceGA)
library(JuliaCall)
library(raster)
library(terra)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(rnaturalearth)
library(rnaturalearthhires)
library(ggplot2)
library(ggspatial)
library(sf)
library(ggtext)
library(scico)

# ---------- #
# Prepare Light Raster ####
# ---------- #

# Read in coordinates
coords <- read.csv("../Coordinates/coordinates.csv")
coords_vect <- vect(coords, geom = c("Lon", "Lat"), crs = "epsg:4326")
plot(coords_vect)

# Coordinates extent
extent <- ext(c(-5.9, -1, 49.9, 53.5))

# Read in nightlight geotiff (from Orly)
nightlight <- rast("data/WGS/spatial_data/light_tiff.tif")
names(nightlight) <- "Nightlight"

# Crop light raster
nightlight <- crop(nightlight, extent)
plot(nightlight)

# UK coastline geopackage
uk_coastline <- vect("data/WGS/spatial_data/UK_borders.gpkg")
plot(uk_coastline)

# Mask sea in light raster
nightlight <- mask(nightlight, uk_coastline)
plot(nightlight)

# Rescale light raster to have values between 1 and 100
r_min <- global(nightlight, "min", na.rm = TRUE)[[1]]
r_max <- global(nightlight, "max", na.rm = TRUE)[[1]]
nightlight_scaled <- ((nightlight - r_min) / (r_max - r_min)) * 99 + 1
plot(nightlight_scaled)
points(coords_vect, col = "white", pch = 21, bg = NA, cex = 1)

# Convert NA (sea cells) to 200
nightlight_scaled[is.na(nightlight_scaled)] <- 200
nightlight_scaled; plot(nightlight_scaled)

# ---------- #
# Prepare Land Cover Rasters ####
# ---------- #

# Read in coordinates
coords <- read.csv("../Coordinates/coordinates_BritishNationalGrid.csv")
coords_vect <- vect(coords, geom = c("Lon", "Lat"), crs = "epsg:27700")
plot(coords_vect)

# Coordinates extent
extent_bng <- ext(c(-5.9, -1, 49.9, 53.5)) |> project("epsg:4326", "epsg:27700")

# Read in 1km Land Cover 2023
land_cover <- rast("data/WGS/spatial_data/gb2023lcm1km_dominant_aggregate.tif")
plot(land_cover)

# Crop land cover to study area
land_cover <- crop(land_cover, extent_bng)
plot(land_cover)

# Read in Excel file with land cover scenario costings
costings <- readxl::read_xlsx("data/WGS/spatial_data/land_cover_scenario_costings.xlsx", n_max = 11)
costings

# Categories:
# 0 = Sea
# 1 = Broadleaf woodland
# 2 = Coniferous woodland
# 3 = Arable
# 4 = Improved grassland
# 5 = Semi-natural grassland
# 6 = Mountain, heath and bog
# 7 = Saltwater
# 8 = Freshwater
# 9 = Coastal
# 10 = Built-up areas and gardens

# Function to assign costings to new land cover rasters
assign_costings <- function (raster, costings_df, scenario) {
  # subset costings data.frame by scenario
  df <- dplyr::select(costings_df, Class_ID, ends_with(as.character(scenario)))

  # create new raster in WGS84 projection to classify new costings
  costings_ras <- project(raster, "EPSG:4326", method = "near")

  # reclassify the raster using the costings lookup table
  costings_ras <- classify(costings_ras, as.matrix(df))
  names(costings_ras) <- ifelse(is.numeric(scenario), paste0("Land_Cover_Scenario", scenario), paste0("Land_Cover_", scenario))

  return(costings_ras)
}

# Create new rasters for costing scenarios
land_cover_S1 <- assign_costings(land_cover, costings, 1)
land_cover_S2 <- assign_costings(land_cover, costings, 2)
land_cover_S3 <- assign_costings(land_cover, costings, 3)
land_cover_Null <- assign_costings(land_cover, costings, "Null")

# Plot new rasters to inspect reclassification
plot(land_cover_S1)
plot(land_cover_S2)
plot(land_cover_S3)
plot(land_cover_Null)

# ---------- #
# Prepare National Forest Inventory Raster ####
# ---------- #

# Preprocessing
# forest <- vect("../../../Downloads/National_Forest_Inventory_GB_2023/National_Forest_Inventory_GB_2023.shp")
# forest <- crop(forest, extent_bng)
# forest <- forest[forest$CATEGORY == "Woodland", ]
# r <- rast(forest, resolution = 1000)
# woodland_raster <- rasterize(forest, r, field = 1, background = NA)
# writeRaster(woodland_raster, filename = "data/WGS/spatial_data/National_Forest_Inventory_GB_2023_study_area_1km.tif", overwrite = TRUE)

# Read in National Forest Inventory 2023 raster
forest <- rast("data/WGS/spatial_data/National_Forest_Inventory_GB_2023_study_area_1km.tif")
plot(forest)

# Compute Euclidean distance to woodland
forest_dist <- terra::distance(forest, method = "geo", unit = "km")
names(forest_dist) <- "Forest_distance"
forest_dist; plot(forest_dist)

# Reproject to WGS84
forest_dist <- project(forest_dist, "EPSG:4326")
plot(forest_dist)

# Mask sea in forest raster
forest_dist <- mask(forest_dist, uk_coastline)
plot(forest_dist)

# Rescale forest raster to have values between 1 and 100
r_min <- global(forest_dist, "min", na.rm = TRUE)[[1]]
r_max <- global(forest_dist, "max", na.rm = TRUE)[[1]]
forest_scaled <- ((forest_dist - r_min) / (r_max - r_min)) * 99 + 1
plot(forest_scaled)

# Convert NA (sea cells) to 200
forest_scaled[is.na(forest_scaled)] <- 200
forest_scaled; plot(forest_scaled)

# ---------- #
# Circuitscape ####
# ---------- #

# Function to automate executation of Circuitscape on rasters
runCircuitscape <- function (JULIA_PATH, OUTPUT_DIR, rast, coordinates) {
    current_map <- ResistanceGA::Run_CS.jl(
      r = raster::raster(rast),
      CS_Point.File = coordinates,
      CurrentMap = TRUE,
      EXPORT.dir = OUTPUT_DIR,
      output = "raster",
      cholmod = TRUE,
      parallel = TRUE,
      cores = 4,
      JULIA_HOME = JULIA_PATH,
      Julia_link = "JuliaCall"
  )
  return(current_map)
}

# List of rasters to run with Circuitscape
cs_rasters <- list(nightlight_scaled, land_cover_S1, land_cover_S2, land_cover_S3, land_cover_Null, forest_scaled)

# Make a directory called Circuitscape in outputs/
# If true (re)run Circuitscape on rasters specified
runCS <- FALSE
JULIA_HOME = "C:/Users/tj311/.julia/juliaup/julia-1.11.5+0.x64.w64.mingw32/bin"
OUTPUT_DIR = "outputs/Circuitscape/"
cs_coords = "../Coordinates/coordinates.txt"

# Run Circuitscape
if (runCS) {
  cs_results <- lapply(cs_rasters, \(x) runCircuitscape(JULIA_HOME, OUTPUT_DIR, x, cs_coords))
}

# ---------- #
# Genetic Distances ####
# ---------- #

# Load dartR
library(dartR)

# Read in SNPs
load("data/WGS/snps_genlight.RData")
snps

# Compute Fst between sites
runFst = FALSE
if (runFst) {
  fst <- dartR::gl.fst.pop(snps, nclusters = 4)
  fst_df <- drop_na(as.data.frame(as.table(fst)))
  names(fst_df) <- c("Site1", "Site2", "Fst")
  write.csv(
    fst_df,
    file = "outputs/WGS_snps_Fst_dataframe.csv",
    row.names = FALSE
  )
}

# Read in Fst data.frame
Fst <- read.csv("outputs/WGS_snps_Fst_dataframe.csv")

# ---------- #
# Euclidean Distances Between Sites ####
# ---------- #

# Read in coordinates in British National Grid and convert to SpatVector object
coords_bng <- "../Coordinates/coordinates_BritishNationalGrid.csv" |> 
  read.csv() |> 
  vect(geom = c("Lon", "Lat"), crs = "epsg:27700")

# Compute Euclidean pairwise distances
euclidean_dist <- terra::distance(coords_bng, unit = "km", pairs = TRUE)
colnames(euclidean_dist) <- c("Site1","Site2","Distance_km")

# Export as CSV in column format
write.csv(euclidean_dist, "outputs/WGS_euclidean_distances.csv", row.names = FALSE)

# Read in all data and build data.frame for modelling
euclidean_dist <- read.csv("outputs/WGS_euclidean_distances.csv")
fst_df <- read.csv("outputs/WGS_snps_Fst_dataframe.csv")
nightlight_cs_resist <- read.delim("outputs/Circuitscape/Nightlight_resistances_3columns.out", header = FALSE, sep = " ")
forest_dist_cs_resist <- read.delim("outputs/Circuitscape/Forest_distance_resistances_3columns.out", header = FALSE, sep = " ")
lc1_cs_resist <- read.delim("outputs/Circuitscape/Land_Cover_Scenario1_resistances_3columns.out", header = FALSE, sep = " ")
lc2_cs_resist <- read.delim("outputs/Circuitscape/Land_Cover_Scenario2_resistances_3columns.out", header = FALSE, sep = " ")
lc3_cs_resist <- read.delim("outputs/Circuitscape/Land_Cover_Scenario3_resistances_3columns.out", header = FALSE, sep = " ")
lcnull_cs_resist <- read.delim("outputs/Circuitscape/Land_Cover_Null_resistances_3columns.out", header = FALSE, sep = " ")
modelling_data <- data.frame(
  Site1 = euclidean_dist$Site1,
  Site2 = euclidean_dist$Site2,
  euclidean_dist = euclidean_dist[[3]],
  euclidean_dist_log = log(euclidean_dist[[3]]),
  fst = fst_df[[3]],
  nightlight_cs_resist = nightlight_cs_resist[[3]],
  forest_dist_cs_resist = forest_dist_cs_resist[[3]],
  lc1_cs_resist = lc1_cs_resist[[3]],
  lc2_cs_resist = lc2_cs_resist[[3]],
  lc3_cs_resist = lc3_cs_resist[[3]],
  lcnull_cs_resist = lcnull_cs_resist[[3]]
)
head(modelling_data)

# Export modelling data table
write.csv(modelling_data, file = "outputs/WGS_landscape_genetics_modelling_data.csv", row.names = FALSE)

# ---------- #
# Linear Mixed-Effects Models ####
# ---------- #

# Load modelling packages
library(lme4)
library(usdm)
library(ecodist)
library(performance)

# Read in data if needed
modelling_data <- read.csv("outputs/WGS_landscape_genetics_modelling_data.csv")

# Test for IBD using Fst and Euclidean distances
MRM(modelling_data$fst ~ modelling_data$euclidean_dist_log, nperm = 10000, mrank = FALSE)

# MRM model R2 > 0.50 with euclidean distance
# Used a formula to control for the effect of Isolation-By-Distance on Fst
modelling_data$fst_ibd = (log(modelling_data$fst+1)) / (modelling_data$euclidean_dist_log)

# Single regression models on distance matrices
# (mrm_light <- MRM(modelling_data$fst_ibd ~ modelling_data$nightlight_cs_resist, nperm = 1000, mrank = FALSE))
# (mrm_forest <- MRM(modelling_data$fst_ibd ~ modelling_data$forest_dist_cs_resist, nperm = 1000, mrank = FALSE))
# (mrm_lc1 <- MRM(modelling_data$fst_ibd ~ modelling_data$lc1_cs_resist, nperm = 1000, mrank = FALSE))
# (mrm_lc2 <- MRM(modelling_data$fst_ibd ~ modelling_data$lc2_cs_resist, nperm = 1000, mrank = FALSE))
# (mrm_lc3 <- MRM(modelling_data$fst_ibd ~ modelling_data$lc3_cs_resist, nperm = 1000, mrank = FALSE))
# (mrm_lcnull <- MRM(modelling_data$fst_ibd ~ modelling_data$lcnull_cs_resist, nperm = 1000, mrank = FALSE))

# Test for multicollinearity after removing non-significant variables (VIF < 5)
modelling_data_sig <- dplyr::select(modelling_data, nightlight_cs_resist:lcnull_cs_resist)
usdm::vif(modelling_data_sig)

# Remove scenario 3 and re-test
modelling_data_sig <- dplyr::select(modelling_data_sig, -lc3_cs_resist)
usdm::vif(modelling_data_sig)

# Remove scenario 2 and re-test
modelling_data_sig <- dplyr::select(modelling_data_sig, -lc2_cs_resist)
usdm::vif(modelling_data_sig)

# Remove scenario null and re-test
modelling_data_sig <- dplyr::select(modelling_data_sig, -lcnull_cs_resist)
usdm::vif(modelling_data_sig)
head(modelling_data_sig)

# Add information to modelling_data_sig data.frame
modelling_data_sig <- cbind(
  modelling_data[,1:2],
  modelling_data_sig,
  fst = modelling_data$fst,
  # MRM model R2 > 0.50 with euclidean distance
  # Used a formula to control for the effect of Isolation-By-Distance
  fst_ibd = (log(modelling_data$fst+1)) / (modelling_data$euclidean_dist_log)
)
head(modelling_data_sig)

# Maximum Likelihood Population Effects (MLPE)
# Models pairwise genetic distances as a function of landscape variables
# Accounts for the non-independence of pairwise data `(1|Site1)`

# Run MLPE with no REML
forest <- lmer(fst_ibd ~ forest_dist_cs_resist + (1|Site1), data = modelling_data_sig, REML = FALSE)
forest_light <- lmer(fst_ibd ~ forest_dist_cs_resist + nightlight_cs_resist + (1|Site1), data = modelling_data_sig, REML = FALSE)
forest_light_lc1 <- lmer(fst_ibd ~ forest_dist_cs_resist + nightlight_cs_resist + lc1_cs_resist + (1|Site1), data = modelling_data_sig, REML = FALSE)
forest_lc1 <- lmer(fst_ibd ~ forest_dist_cs_resist + lc1_cs_resist + (1|Site1), data = modelling_data_sig, REML = FALSE)
light <- lmer(fst_ibd ~ nightlight_cs_resist + (1|Site1), data = modelling_data_sig, REML = FALSE)
light_lc1 <- lmer(fst_ibd ~ nightlight_cs_resist + lc1_cs_resist + (1|Site1), data = modelling_data_sig, REML = FALSE)
lc1 <- lmer(fst_ibd ~ lc1_cs_resist + (1|Site1), data = modelling_data_sig, REML = FALSE)

# Evaluate models with AICc and BIC
performance::compare_performance(
  forest, forest_light, forest_light_lc1, light, light_lc1, lc1, forest_lc1,
  estimator = "ML",
  rank = TRUE
)

# ---------- #
# Visualisation ####
# ---------- #

# Plot Fst versus Distance to Forest Resistance
(plt_Fst_Forest <- ggplot(modelling_data_sig, aes(x = forest_dist_cs_resist, y = fst_ibd))+
  geom_point()+
  scale_y_continuous(
    limits = c(0, 0.03),
    expand = c(0,0),
  )+
  labs(
    x = "Distance to Woodland Resistance",
    y = "<span style='font-size:10pt;'>Genetic Differentiation</span><br><span style='font-size:8pt'>log(<em>F</em><sub>st</sub> + 1) / log(Euclidean Distance)</span>",
    tag = "A"
   )+
  theme(
    panel.background = element_rect(fill = "#F5F5F5"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
    axis.title.y = element_markdown(),
    axis.title.x = element_text(size = 10, vjust = 0),
    axis.text = element_text(size = 6),
  )
)

# Plot Fst versus Nightlight Resistance
(plt_Fst_Nightlight <- ggplot(modelling_data_sig, aes(x = nightlight_cs_resist, y = fst_ibd))+
  geom_point()+
  scale_y_continuous(
    limits = c(0, 0.03),
    expand = c(0,0),
  )+
  labs(
    x = "ALAN Resistance",
    y = "<span style='font-size:10pt; font-weight: bolder;'>Genetic Differentiation</span><br><span style='font-size:8pt'>log(<em>F</em><sub>st</sub> + 1) / log(Euclidean Distance)</span>",
    tag = "B"
   )+
  theme(
    panel.background = element_rect(fill = "#F5F5F5"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
    axis.title.y = element_markdown(),
    axis.title.x = element_text(size = 10, vjust = 0),
    axis.text = element_text(size = 6),
  )
)

# Read in coordinates
coords <- read.csv("../Coordinates/coordinates.csv")
coords_vect <- vect(coords, geom = c("Lon", "Lat"), crs = "epsg:4326")

# Read in Circuitscape cumulative current rasters
current_forest <- rast("outputs/Circuitscape/Forest_distance_cum_curmap.asc")
current_light <- rast("outputs/Circuitscape/Nightlight_cum_curmap.asc")

# Remove sea from rasters
current_forest <- mask(current_forest, uk_coastline)
current_light <- mask(current_light, uk_coastline)

# Basic plot()
plot(current_forest)
plot(current_light)

# Get the maximum current value outside of site coordinate cells
forest_cells_to_exclude <- cells(current_forest, coords_vect)
forest_cells_to_exclude <- forest_cells_to_exclude[!is.na(forest_cells_to_exclude)]
max_val_forest <- max(values(current_forest)[-forest_cells_to_exclude], na.rm = TRUE)

# Get the third quartile value 
light_cells_to_exclude <- cells(current_light, coords_vect)
light_cells_to_exclude <- light_cells_to_exclude[!is.na(light_cells_to_exclude)]
values(current_light) |> summary()
# max_val_light <- max(values(current_light)[-light_cells_to_exclude], na.rm = TRUE)
max_val_light <- 1

# Rescale site cells to maximum value so that current patterns are easier to interpret
current_forest[current_forest > max_val_forest] <- max_val_forest
current_light[current_light > max_val_light] <- max_val_light

# Check rescaling using basic plot()
plot(current_forest)
plot(current_light)

# Reproject spatial data
current_forest_bng <- project(current_forest, "EPSG:27700")
current_light_bng <- project(current_light, "EPSG:27700")
uk_coastline_bng <- project(uk_coastline, "EPSG:27700")

# Colour palette
library(RColorBrewer)
library(viridis)
col_palette <- "Oranges"
max_colours <- brewer.pal.info[col_palette, ]$maxcolors
# max_colours <- 8
# viridis::plasma(10)
fill_colours <- viridis::viridis(5)

# Extent
extent <- project(ext(-5.8, -1.1, 49.9, 53.1), "EPSG:4326", "EPSG:27700")

# Custom theme
map_theme <- theme(
  panel.background = element_rect(fill = "#F5F5F5"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
  axis.title.x = element_text(vjust = 0, size = 10),
  axis.title.y = element_text(size = 10),
  axis.text = element_text(size = 6),
  # legend.title = element_text(size = 10, vjust = 3),
  legend.title = ggtext::element_markdown(size = 10, margin = margin(b = 10)),
  legend.text = element_text(size = 9),  
)

# Function to plot current map
plot_current_map <- function (var, rast, tag, max_val) {
  ggplot()+
    ggspatial::layer_spatial(rast)+
    ggspatial::layer_spatial(data = uk_coastline_bng, fill = NA)+
    scico::scale_fill_scico(
      name = "<span>Movement<br>Density<br>Potential<span>",
      palette = "vik",
      na.value = NA,
      limits = c(0, max_val),
      breaks = c(0, max_val),
      labels = c("Low", "High")
    )+
    labs(
      x = "Longitude",
      y = "Latitude",
    )+
    annotate(
      "text",
      x = extent$xmax, y = extent$ymax,
      hjust = 1.05, vjust = 1.60, fontface = "bold",
      label = var, colour = "white", size = 3.5
    )+
    coord_sf(
      xlim = c(extent$xmin, extent$xmax),
      ylim = c(extent$ymin, extent$ymax),
      expand = FALSE
    )+
    labs(
      tag = tag
    )+
    map_theme
}

# Plot distance to forest current map
(plt_map_forest <- plot_current_map("Distance to Woodland", current_forest_bng, "C", max_val_forest))

# Plot nightlight current map
(plt_map_light <- plot_current_map("ALAN ", current_light_bng, "D", max_val_light))

# Prepare figure using patchwork
library(patchwork)

# Figure composition
layout <- "
AB
CD
"
wrap_plots(
  A = plt_Fst_Forest,
  B = plt_Fst_Nightlight+ theme(axis.title.y = element_blank()),
  C = plt_map_forest+ theme(
    legend.position = "none",
    axis.title.y = element_text(vjust = 5), 
  ),
  D = plt_map_light+ theme(
    plot.margin = margin(l = 30, unit = "pt"),
    legend.position = "right", legend.margin = margin(l = 10, unit = "pt"),
    axis.title.y = element_blank(), 
  ),
  design = layout
)
ggsave("figures/Figure4.png", width = 10, height = 8, units = "in", dpi = 600)
ggsave("figures/Figure4.pdf", width = 10, height = 8, units = "in")
