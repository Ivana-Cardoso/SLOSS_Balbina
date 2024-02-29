# SLOSS cube hypothesis test
# in Balbina Hydroelectric Dam landscapes
#
# Author: Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: February 28, 2024

set.seed(13) # L

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")

# Load packages
library(raster)
library(sf)
library(terra)
library(landscapemetrics)
library(ggplot2)
library(ggpubr)

#### SELECTING LANDSCAPE SITES ####
# I downloaded the land-use cover raster for Brazil in the year 2022 using the MAPBIOMAS Toolkit (https://brasil.mapbiomas.org/colecoes-mapbiomas/) on the Google Engine. In the Google Engine itself, I selected the area of interest by using a shapefile as a mask. This shapefile is available on my GitHub (https://github.com/Ivana-Cardoso/SLOSS_Balbina) as 'Balbina.shp' if you wish to crop it as I did. Alternatively, the raster file is also available on my GitHub as 'Balbina_mapbiomas_collection8_2022.tif'.

# Import data
Balbina <- raster("Balbina_mapbiomas_collection8_2022.tif")
new_crs <- CRS("+proj=utm +south +zone=21 +datum=WGS84 +units=m +no_defs")
Balbina <- projectRaster(Balbina, crs=new_crs)
values(Balbina)  = as.integer(values(Balbina))
check_landscape(Balbina)
plot(Balbina)

# Below, I will (1) create an object only with forest pixels to calculate forest cover and the number of forest patches (habitat), and (2) create a new raster with all land-use types except forest (non-habitat)

# Reference code for land-use type pixel values:  https://brasil.mapbiomas.org/wp-content/uploads/sites/4/2023/08/Legenda-Colecao-8-LEGEND-CODE-1.pdf

# Select only forest pixels 
forest_formation <- Balbina == 3
forest_formation <- projectRaster(forest_formation, crs=new_crs)
values(forest_formation)  = as.integer(values(forest_formation))
check_landscape(forest_formation)
plot(forest_formation)

# Select all land-use types except forest
non_habitat <- Balbina
non_habitat[!(Balbina != 3)] <- NA
non_habitat <- projectRaster(non_habitat, crs=new_crs)
values(non_habitat) <- as.integer(values(non_habitat))
check_landscape(non_habitat)
plot(non_habitat)

# Create a grid of points with a 4km interval
extent_grid <- as(extent(forest_formation), 'SpatialPolygons')
proj4string(extent_grid) <- crs(new_crs)
points <- sp::makegrid(extent_grid, cellsize = 4000)  # Evenly spaced points (4km)
plot(forest_formation)
points(points, pch = ".")

# Calculate the number of patches in a 2 km-radius landscape
# Reference: https://r-spatialecology.github.io/landscapemetrics/articles/guide_sample_lsm.html
grid <- st_as_sf(points, coords = c("x1", "x2"), crs = new_crs)
number_patches <- landscapemetrics::sample_lsm(forest_formation, y = grid, 
                                               size = 2000,
                                               shape = "square",
                                               what = "lsm_c_np",
                                               directions = 8)
number_patches <- subset(number_patches, number_patches$class == 1)  # Select only the presence of patches

# Calculate forest cover in a 2 km-radius landscape
forest_cover <- landscapemetrics::sample_lsm(forest_formation, y = grid, size = 2000,
                                             shape = "square",
                                             what = "lsm_c_pland",
                                             directions = 8)
forest_cover <- subset(forest_cover, forest_cover$class == 1)  # Select only the presence of forest

# Plot number of patches X forest cover
number_patches$plot_id == forest_cover$plot_id
metrics <- as.data.frame(cbind(number_patches$plot_id, number_patches$value,
                               forest_cover$value))
colnames(metrics) <- c("id", "number_patches", "forest_cover")

NP_FC_plot <- 
  ggplot(data = metrics,
         mapping = aes(x = number_patches, y = forest_cover)) +
  labs(x = "Number of forest patches",
       y = "Forest cover (%)") +
  geom_point() +
  geom_hline(yintercept = c(25, 50))+
  theme_pubr(base_size = 20)

NP_FC_plot 

ggsave("Number_patches_vs_Forest_cover.png", plot=NP_FC_plot,
       width = 200, height = 200, units = "mm",
       dpi = 600)

# Select landscapes with 25 to 50% forest cover because, within these values, forest cover is visually not correlated with number of forest fragments
metrics = subset(metrics, metrics$forest_cover <= 50 & metrics$forest_cover >= 25)
cor(metrics$number_patches, metrics$forest_cover) # not correlated (r=-0.14)

# Within this selection (from 25% to 50% forest cover), I searched the NP_FC_plot graph for points at the extremes of the number of patches gradient to make clear the characterization of few large patches and several small patches. I checked the distribution of each point on the map and sought to make a broad selection that covered a significant portion of the study area. I discarded points that I deemed inaccessible due to safety concerns. Below are all selected and discarded landscapes identified by their IDs.
selected_landscapes <- subset(metrics, 
                              metrics$id %in% c(145, 583, 503, 112, 428,
                                                259, 468, 229, 422, 381,
                                                188, 228, 28, 226, 186,
                                                545, 585, 146, 586, 546, 
                                                76, 299, 588, 387, 427, 
                                                150, 384, 542, 549, 338,
                                                301, 584))

discarded_landscapes <- subset(metrics, metrics$id %in% c(690, 727, 646, 654,
                                                          789, 728, 489, 331,
                                                          749, 693, 687, 614,
                                                          452, 624, 604))

NP_FC_plot +
  geom_point(selected_landscapes, 
             mapping = aes(x=number_patches, y=forest_cover), col = "blue")+
  geom_point(discarded_landscapes, 
             mapping = aes(x=number_patches, y=forest_cover), col = "red")

cor(selected_landscapes$number_patches, selected_landscapes$forest_cover) # checking correlation of number of patches and forest cover in selected landscapes. Not correlated (r=-0.058)

# Visualize in map
points$id <- seq(1,nrow(points), 1)
landscapes <- st_buffer(grid, endCapStyle = "SQUARE", dist = 2000)
landscapes$id <- points$id
landscapes <- landscapes[landscapes$id %in% selected_landscapes$id,]
landscapes <- merge(landscapes, metrics, by = "id")
plot(forest_formation)
plot(landscapes$geometry, add = TRUE, col = "blue")

# I need 10 landscapes with few large patches and 10 landscapes with many small patches. Therefore, I selected 10 of each, avoiding neighboring landscapes, and maintaining SS above 40 fragments and SL below 20.
landscapes <- landscapes[-c(5, 8, 9, 11, 13, 20, 21, 24, 26, 27, 28, 30),]
plot(forest_formation)
plot(landscapes$geometry, add=TRUE, col="blue")

cor(landscapes$number_patches, landscapes$forest_cover) # Not correlated (r=-0.019)

# Mean and sd number of forest fragments and forest cover for our 'Several Small' and 'Single Large' landscapes
SS =  landscapes[landscapes$number_patches >= 40,]
SL = landscapes[landscapes$number_patches <= 20,]

mean(SS$number_patches) # 45.9
sd(SS$number_patches) # 7.17
mean(SS$forest_cover) # 35.86
sd(SS$forest_cover) # 7.65

mean(SL$number_patches) # 14.9
sd(SL$number_patches) # 3.28
mean(SL$forest_cover)# 35.11
sd(SL$forest_cover) # 7.75


# Calculate the coverage of other land-uses composing the selected landscapes besides the forest formation.
non_habitat_cover = landscapemetrics::sample_lsm(non_habitat, 
                                                 y = landscapes$geometry,
                                                 plot_id = landscapes$id,
                                                 shape = "square",
                                                 what = "lsm_c_pland",
                                                 directions = 8)

# Export files
dir.create("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/Exported_files")
setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/Exported_files")
writeRaster(forest_formation, "forest_formation", format="GTiff", overwrite=T)
writeRaster(non_habitat, "non_habitat", format="GTiff", overwrite=T)
st_write(landscapes, "landscapes.shp", append = F)
write.csv(non_habitat_cover, "non_habitat_cover.csv")


#### SELECTING SAMPLING POINTS FOR DEPLOYMENT OF AUDIO RECORDERS ####
# I modified Pavel's function (with his permission) to randomly select only 10 points in each landscape while respecting a minimum distance of 200 meters. The original function is available at: https://github.com/pdodonov/sampling
setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")
source("siteSelection_v0.8.R")

# Each audio recorder covers an area of 50m, so, to avoid that audio recorder points fall on the edge of the landscape and on the edge of fragments, I will (1) create a buffer 100 m inside each landscape and (2) create a buffer 100 m inside each forest fragment, so the area of the audio recorder will be inside the fragment and not in its edges.
buffer_landscapes <- sf::st_buffer(landscapes, dist = -100)

# Using QGIS v.3.34 I transformed forest_cover raster into polygons, selected only the forest polygons (raster value = 1) and created a buffer 50 m inside all fragments. I will import it:
buffer_forest <- read_sf("Forest_buffer.shp")
buffer_forest <- st_set_crs(buffer_forest, new_crs)

recorder_points <- data.frame()

for (i in 1:20) {
  set.seed(13)
  
  # Extract the values from the extent of landscape i
  landscape <- as(extent(buffer_landscapes[i, 4]), 'SpatialPolygons')
  proj4string(landscape) <- crs(new_crs)
  
  # Create points equally spaced at 100 meters within the extent of landscape i  
  recorders.i <- sp::makegrid(landscape, cellsize = 100)
  recorders.i <- st_as_sf(recorders.i, coords = c("x1", "x2"), crs = new_crs)
  
  # Extracts the raster values at the point locations and selects only the points in the forest (value = 1)
  recorders.i <- recorders.i[buffer_forest,]
  recorders.i <- st_set_crs(recorders.i, new_crs)
  
  # Converts spatial points object into a data.frame and identifies each one.  
  recorders.i <- st_coordinates(recorders.i)
  recorders.i <- as.data.frame(recorders.i)
  recorders.i$id <- 1:nrow(recorders.i)
  
  # Select random points at least 200 meters apart  
  recorders <- select.site(x = recorders.i, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)
  
  # Converts the results into a data.frame, identifies, and adds the coordinates
  recorders <- as.data.frame(recorders)
  colnames(recorders) <- "id"
  recorders <- merge(recorders, recorders.i, by = "id")
  
  recorder_points <- rbind(recorder_points, recorders)
}

plot(forest_formation)
points(recorder_points$X, recorder_points$Y, pch = 16, col = "blue")

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/Exported_files")
# write.csv(recorder_points, "recorder_points.csv")

# Out of curiosity: fragments size mean/max/min/sd
setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/Patches")
SS_area <- read_sf("Patches_SS.shp")
SS_area <- SS_area$hectares
min(SS_area)
max(SS_area)
mean(SS_area)
sd(SS_area)
hist(SS_area)

SL_area <- read_sf("Patches_SL.shp")
SL_area <- SL_area$hectares
SL_area <- SL_area[-c(19:21, 23, 27, 31)] # removing continuous forest
min(SL_area)
max(SL_area)
mean(SL_area)
sd(SL_area)
hist(SL_area)

