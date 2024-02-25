# SLOSS cube hypothesis test
# in Balbina Hydroelectric Dam landscapes
#
# Author: Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: February 24, 2024

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")

# Load packages
library(raster)
library(sf)
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



