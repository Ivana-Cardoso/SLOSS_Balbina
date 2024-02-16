# SLOSS cube hypothesis test
# in Balbina Hidroeletric Dam landscapes
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: February 16, 2024

# Set working directory
setwd("C:/Users/Ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")

# Load packages
library(raster)
library(terra)
library(sf)
library(landscapemetrics)

# Import data
MAPBIOMAS <- raster("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_8/lclu/coverage/brasil_coverage_2022.tif")
proj4string(MAPBIOMAS) <- CRS("+init=epsg:4326")
plot(MAPBIOMAS)


# Generate data for Balbina
Balbina <- as(extent(-60.6,-59.2,-2.1, -0.9), 'SpatialPolygons') # creating crop area that corresponds to Balbina Hydroeletric Dam
proj4string(Balbina) = crs("+init=epsg:4326")

extent(MAPBIOMAS)
extent(Balbina)
plot(Balbina, add=TRUE)

Balbina_raster <- raster::crop(MAPBIOMAS, Balbina) # cropping mapbiomas by Balbina
new_crs <- CRS("+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs")
Balbina_raster <- projectRaster(Balbina_raster, crs=new_crs)
landscapemetrics::check_landscape(Balbina_raster)
plot(Balbina_raster)


# Select only pixels of forest
Forest_formation <- Balbina_raster == 3 # selecting only pixel values that correspond to forest. More information in: https://brasil.mapbiomas.org/wp-content/uploads/sites/4/2023/08/Legenda-Colecao-8-LEGEND-CODE-1.pdf
Forest_formation <- projectRaster(Forest_formation, crs=new_crs)
Forest_formation = round(Forest_formation)
landscapemetrics::check_landscape(Forest_formation)
plot(Forest_formation)


# Create grid of points with 4km interval
extent_grid <- as(extent(Forest_formation), 'SpatialPolygons')
grid <- sp::makegrid(extent_grid, cellsize = 4000) # evenly spaced points (4km)
points(grid, pch=".")


