# SLOSS cube hypothesis test
# in Balbina Hidroeletric Dam landscapes
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: February 19, 2024

# Set working directory
setwd("C:/Users/Ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")

# Load packages
library(raster)
library(terra)
library(sf)
library(landscapemetrics)
library(ggplot2)
library(ggpubr)

#### SELECTING LANDSCAPE SITES ####

# Import data
MAPBIOMAS <- raster("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_8/lclu/coverage/brasil_coverage_2022.tif")
proj4string(MAPBIOMAS) <- CRS("+init=epsg:4326")
plot(MAPBIOMAS)


# Generate data for Balbina
Balbina <- as(extent(-60.6,-59.2,-2.1, -0.9), 'SpatialPolygons') # creating a cropping area corresponding to Balbina Hydroelectric Dam
proj4string(Balbina) = crs("+init=epsg:4326")

extent(MAPBIOMAS)
extent(Balbina)
plot(Balbina, add=TRUE)

Balbina_raster <- raster::crop(MAPBIOMAS, Balbina)
new_crs <- CRS("+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs")
Balbina_raster <- projectRaster(Balbina_raster, crs=new_crs)
landscapemetrics::check_landscape(Balbina_raster)
plot(Balbina_raster)


# Select only forest pixels
Forest_formation <- Balbina_raster == 3 # select only pixel values corresponding to forest. Reference: https://brasil.mapbiomas.org/wp-content/uploads/sites/4/2023/08/Legenda-Colecao-8-LEGEND-CODE-1.pdf
Forest_formation <- projectRaster(Forest_formation, crs=new_crs)
Forest_formation = round(Forest_formation)
landscapemetrics::check_landscape(Forest_formation)
plot(Forest_formation)


# Create grid of points with 4km interval
extent_grid <- as(extent(Forest_formation), 'SpatialPolygons')
grid <- sp::makegrid(extent_grid, cellsize = 4000) # evenly spaced points (4km)
points(grid, pch=".")


# Calculate number of forest patches in 2 km-radius landscape
# https://r-spatialecology.github.io/landscapemetrics/articles/guide_sample_lsm.html
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)
number_patches = landscapemetrics::sample_lsm(Forest_formation, y = grid, size=2000,
                                              what = "lsm_c_np",
                                              directions = 8)
number_patches = subset(number_patches, number_patches$class == 1)


# Calculate forest cover in 2 km-radius landscape
forest_cover = landscapemetrics::sample_lsm(Forest_formation, y = grid, size=2000,
                                            what = "lsm_c_pland",
                                            directions = 8)
forest_cover = subset(forest_cover, forest_cover$class == 1)


# Plot number of patches X forest cover
number_patches$plot_id == forest_cover$plot_id
metrics <- as.data.frame(cbind(number_patches$plot_id, number_patches$value,
                               forest_cover$value))
colnames(metrics) = c("id", "number_patches", "forest_cover")

NP_FC_plot = 
  ggplot(data = metrics,
       mapping = aes(x = number_patches, y = forest_cover)) +
  labs(x = "Number of forest patches",
       y = "Forest cover (%)") +
  geom_point() +
  theme_pubr(base_size = 20)
ggsave("Number_patches_vs_Forest_cover.png", plot=NP_FC_plot,
       width = 200, height = 200, units = "mm",
       dpi = 600)


