# SLOSS cube hypothesis test
# in Balbina Hidroeletric Dam landscapes
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: February 13, 2024

# Set working directory
setwd("C:/Users/Ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")

# Load packages
library(raster)
library(terra)


# Import data
MAPBIOMAS <- raster("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_8/lclu/coverage/brasil_coverage_2022.tif")
crs(MAPBIOMAS) <- CRS("+init=epsg:4326")
plot(MAPBIOMAS)

# Handling data
Balbina <- as(extent(-61.0847206680000028,-58.3930087510000035,-2.4154799670000000, -0.4072063180000000), 'SpatialPolygons')
crs(Balbina) = crs(MAPBIOMAS)

extent(MAPBIOMAS)
extent(Balbina)

plot(Balbina, add=TRUE)

Balbina_raster <- raster::crop(MAPBIOMAS, Balbina)
plot(Balbina_raster)


