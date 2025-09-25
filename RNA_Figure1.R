# ========== #
#
# Rhipposideros WGS Analysis 2025
#
# Figure 1
#
# ========== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(sf)
library(ggplot2)
library(terra)
library(dplyr)
library(stringr)
library(mapmixture)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthhires)
library(jpeg)
library(grid)
library(patchwork)

# ---------- #
# Sampling Map
# ---------- #

# CRS
CRS = 4326

# Boundary
bbox <- mapmixture::transform_bbox(c(xmin = -6.7, xmax = -1.1, ymax = 49.9, ymin =  53.5), CRS)

# UK coastline geopackage
uk_coastline <- vect("data/WGS/spatial_data/UK_borders.gpkg")
# plot(uk_coastline)

# Read in light raster
light <- rast("data/WGS/spatial_data/light_tiff.tif") |> 
  terra::mask(x = _, uk_coastline) |> 
  # terra::project(x = _, y = str_c("EPSG:", CRS)) |> 
  terra::crop(x = _, y = bbox) 
# plot(light)

# Cheddar
cheddar <- data.frame(
  Site = "CHE",
  County = "Somerset",
  Country = "England",
  Lon = -2.78,
  Lat = 51.28
)

# Read in WGS site coordinates
coords <- read.csv("../Coordinates/coordinates.csv")
coords <- rbind(coords, cheddar)
coords$Sequencing <- ifelse(coords$Site == "CHE", "RNA-seq", "WGS")
coords_sf <- st_as_sf(coords, coords = c("Lon","Lat"), crs = CRS)
# plot(coords_sf)

# World polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

# Orthographic projection (centered at 0° lon, 0° lat)
ortho_crs <- "+proj=ortho +lat_0=40 +lon_0=10 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Convert bbox to polygon
bbox_poly <- st_sf(st_as_sfc(st_bbox(bbox)), crs = CRS)

# Plot world map
(world_map <- ggplot() +
  geom_sf(data = world, fill = "grey70", colour = "black", linewidth = 0.1) +
  geom_sf(data = bbox_poly, fill = NA, colour = "#ff0000", linewidth = 0.5) +
  coord_sf(
    crs = ortho_crs,
    expand = FALSE,
    xlim = c(-2.5e6, 1.5e6),
    ylim = c(-0.8e6, 3.5e6),
  )+
  theme_minimal()+
  theme(
    axis.text = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(linewidth = 0.2)
  )
)

# Custom theme
custom_theme <- theme(
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "#f0f0f0", linewidth = 0.25),
  axis.title.x = element_text(vjust = -0.5, size = 8),
  axis.title.y = element_text(vjust = 3, size = 8),
  axis.text = element_text(size = 6),
  legend.title = element_text(vjust = 2, size = 10),
  legend.key.width = unit(0.40, "cm"),
  legend.spacing.y = unit(1, "cm")
)

# Plot map
(sampling_map <- ggplot()+
  layer_spatial(data = light)+
  geom_sf(
    data = coords_sf,
    aes(shape = Sequencing),
    size = 2.3, fill = "#ff0000", colour = "white", alpha = 0.90
  )+
  scale_shape_manual(
    name = "Sequencing",
    values = c(24,21)
  )+
  scale_fill_viridis_c(
    name = "ALAN",  
    na.value = NA,
  )+
  # annotation_north_arrow(
  #   data = coords_sf, location = "br", which_north = "true",
  #   height = unit(0.6, "cm"), width = unit(0.6, "cm"),
  #   pad_y = unit(0.8, "cm"),
  #   style = north_arrow_orienteering(text_size = 5)
  # )+
  annotation_scale(
    data = coords_sf, location = "br",
    width_hint = 0.2, bar_cols = c("black","white"),
    # height = unit(0.5, "cm"),
    text_cex = 0.5
  )+
  coord_sf(expand = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  custom_theme+
  annotation_custom(
    grob = ggplotGrob(world_map),
    xmin = -10.5, xmax = Inf,
    ymin = 52.35, ymax = Inf
  )
)

# ---------- #
# Bat Image
# ---------- #

# Import jpeg and convert to ggplot object
bat_jpeg <- readJPEG("data/RNA/LHS_flight_Daniel_Whitby.jpeg") |>
  grid::rasterGrob(image = _, interpolate = TRUE)

# Plot image
bat_img <- ggplot()+
  annotation_custom(bat_jpeg, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_void()

# ---------- #
# Fieldwork Image
# ---------- #

# Import jpeg and convert to ggplot object
fieldwork_jpeg <- readJPEG("data/RNA/Fieldwork1_light.jpg") |>
  grid::rasterGrob(image = _, interpolate = TRUE)

# Plot image
fieldwork_img <- ggplot()+
  annotation_custom(fieldwork_jpeg, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_void()+
  annotation_custom(
    bat_jpeg,
    xmin = 0.055, xmax = 0.455, 
    ymin = 0.85, ymax = 1.04
  )

# ---------- #
# Compose Figure
# ---------- #

# Layout design
layout <- "
  AB
"

# Plot layout
plt_list = list(
  sampling_map+ labs(tag = "A"),
  fieldwork_img+ labs(tag = "B")
)
(Figure1 <- wrap_plots(plt_list, design = layout))

# Export
ggsave(plot = Figure1, filename = "figures/Figure1.png", width = 10, height = 5, dpi = 600)
