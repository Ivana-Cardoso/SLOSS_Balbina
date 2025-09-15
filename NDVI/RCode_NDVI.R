# SLOSS cube hypothesis test
# in Balbina Hydroelectric Dam landscapes
#
# Calculating NDVI to assess environmental heterogeneity
#
# Author: Ivana Cardoso - ivanawaters@gmail.com
#
# Created: September 15, 2025
# Last modification: September 15, 2025


## Imagery used from landsat
# https://www.usgs.gov/centers/eros/science/usgs-eros-archive-landsat-archives-landsat-8-9-olitirs-collection-2-level-2
# Landsat 8-9 Operational Land Imager (OLI) and Thermal Infrared (TIRS) Collection 2 Level-2 Science Products 30-meter multispectral data.
# NDVI calculation instructions on: https://www.usgs.gov/landsat-missions/landsat-normalized-difference-vegetation-index

library(terra)
library(sf)

# Importing data
rm(list = ls()); gc()
set.seed(13)
new_crs <- CRS("+proj=utm +south +zone=21 +datum=WGS84 +units=m +no_defs")

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/NDVI/Landsat_L2")

Balbina_red_1 <- rast("LC09_L2SP_231061_20241130_20241201_02_T1_SR_B4.TIF")
Balbina_red_1 <- project(Balbina_red_1, new_crs)
Balbina_red_2 <- rast("LC09_L2SP_230061_20240920_20240923_02_T1_SR_B4.TIF")
Balbina_red_2 <- project(Balbina_red_2, new_crs)

Balbina_nir_1 <- rast("LC09_L2SP_231061_20241130_20241201_02_T1_SR_B5.TIF")
Balbina_nir_1 <- project(Balbina_nir_1, new_crs)
Balbina_nir_2 <- rast("LC09_L2SP_230061_20240920_20240923_02_T1_SR_B5.TIF")
Balbina_nir_2 <- project(Balbina_nir_2, new_crs)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/data")
points <- read.csv("ants_env.csv")
points <- st_as_sf(points, coords = c("long_utm", "lat_utm"), crs = new_crs)

# Calculating NDVI
Balbina_ndvi_1 <- (Balbina_nir_1 - Balbina_red_1) / (Balbina_nir_1 + Balbina_red_1)
plot(Balbina_ndvi_1)

NDVI_1 <- extract(Balbina_ndvi_1, points, xy = TRUE)

Balbina_ndvi_2 <- (Balbina_nir_2 - Balbina_red_2) / (Balbina_nir_2 + Balbina_red_2) 
plot(Balbina_ndvi_2)
NDVI_2 <- extract(Balbina_ndvi_2, points)

merge(NDVI_1, NDVI_2, by = "ID")
