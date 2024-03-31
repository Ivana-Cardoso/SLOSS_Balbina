# SLOSS cube hypothesis test
# in Balbina Hydroelectric Dam landscapes
#
# Author: Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: March 28, 2024
#
# Select landscapes

rm(list = ls()); gc(); dev.off()

set.seed(13)

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
forest_cover <- landscapemetrics::sample_lsm(forest_formation, y = grid, 
                                             size = 2000,
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
metrics1 = subset(metrics, metrics$forest_cover <= 50 & metrics$forest_cover >= 25)
cor.test(metrics1$number_patches, metrics1$forest_cover) # not correlated (r=-0.14)

# Within this selection (from 25% to 50% forest cover), I searched the NP_FC_plot graph for points at the extremes of the number of patches gradient to make clear the characterization of few large patches and several small patches. I checked the distribution of each point on the map and sought to make a broad selection that covered a significant portion of the study area. I discarded points that I deemed inaccessible due to safety concerns. Below are all selected and discarded landscapes identified by their IDs.
# I selected 10 FEW LARGE landscapes, 10 SEVERAL SMALL, 10 landscapes in the middle of the gradient
selected_landscapes <- subset(metrics, 
                              metrics$id %in% c(76, 259, 301, 338, 387, # SS
                                                468, 545, 584, 586, 588, # SS
                                                299, 427, 428, 583, 585, # SS
                                                546, # SS
                                                28, 112, 145, 150, 186, # FL
                                                228, 381, 384, 422, 503, # FL
                                                188, 146, 229, 226, 542, # FL
                                                549, # FL
                                                66, 74, 115, 143, 189, 263, # Middle
                                                269, 307, 340, 506 # Middle
                              ))

discarded_landscapes <- subset(metrics, metrics$id %in% c(690, 727, 646, 654,
                                                          789, 728, 489, 331,
                                                          749, 693, 687, 614,
                                                          452, 624, 604))

# Control landscapes 
# I chose continuous forest landscapes close to where we had already accessed in previous studies
control_landscapes <- subset(metrics, metrics$id %in% c(25, 78, 194, 315, 500))

NP_FC_plot +
  geom_point(selected_landscapes, 
             mapping = aes(x=number_patches, y=forest_cover), col = "blue")+
  geom_point(discarded_landscapes, 
             mapping = aes(x=number_patches, y=forest_cover), col = "red")+
  geom_point(control_landscapes, 
             mapping = aes(x=number_patches, y=forest_cover), col = "green")


cor.test(selected_landscapes$number_patches, selected_landscapes$forest_cover) # checking correlation of number of patches and forest cover in selected landscapes. Not correlated (r=-0.09)

# Visualize in map
points$id <- seq(1,nrow(points), 1)
landscapes <- st_buffer(grid, endCapStyle = "SQUARE", dist = 2000)
landscapes$id <- points$id
landscapes <- landscapes[landscapes$id %in% selected_landscapes$id,]
landscapes <- merge(landscapes, metrics, by = "id")

landscapes_control <- st_buffer(grid, endCapStyle = "SQUARE", dist = 2000)
landscapes_control$id <- points$id
landscapes_control <- landscapes_control[landscapes_control$id %in% control_landscapes$id,]
landscapes_control <- merge(landscapes_control, metrics, by = "id")

# I need 10 landscapes with few large patches, 10 landscapes with many small patches and 10 landscapes in the middle of the gradient (>20 e <40 forest patches). Therefore, I selected them, avoiding neighboring landscapes, and maintaining SS above 40 fragments and SL below 20.
landscapes <- landscapes[-c(1, 9, 12, 14, 16, 20, 29, 30, 34, 36, 38, 40),]
plot(forest_formation)
plot(landscapes$geometry, add=TRUE, col="blue")
plot(landscapes_control$geometry, add = TRUE, col = "green")

cor.test(landscapes$number_patches, landscapes$forest_cover) # Not correlated (r=-0.168)

# Mean and sd number of forest fragments and forest cover for our 'Several Small' and 'Single Large' landscapes
SS =  landscapes[landscapes$number_patches >= 40,]
SL = landscapes[landscapes$number_patches <= 20,]
MS = landscapes[landscapes$number_patches > 20 &  landscapes$number_patches < 40,]

mean(SS$number_patches) # 45.9
sd(SS$number_patches) # 7.17
mean(SS$forest_cover) # 35.86
sd(SS$forest_cover) # 7.65

mean(SL$number_patches) # 15
sd(SL$number_patches) # 3.23
mean(SL$forest_cover)# 36.52
sd(SL$forest_cover) # 8.14

mean(MS$number_patches) # 30.5
sd(MS$number_patches) # 5.06
mean(MS$forest_cover)# 36.83
sd(MS$forest_cover) # 6.87

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
st_write(landscapes_control, "control_landscapes.shp", append = F)
write.csv(non_habitat_cover, "non_habitat_cover.csv")

#### CALCULATING MEAN PATCH SIZE AND EVENNESS IN LANDSCAPES ####
patch_area = landscapemetrics::sample_lsm(forest_formation,
                                          y = landscapes$geometry,
                                          plot_id = landscapes$id,
                                          shape = "square",
                                          what = "lsm_p_area",
                                          directions = 8)
patch_area <- subset(patch_area, patch_area$class == 1)
patch_area <- as.data.frame(cbind(patch_area$plot_id,
                                  patch_area$id,
                                  patch_area$value))
colnames(patch_area) <- c("landscape_id", "patch_id", "area")
patch_area$landscape_id = as.factor(patch_area$landscape_id)

# Mean patch size
mean_patch_size <- aggregate(area ~ landscape_id, data = patch_area, FUN = mean)
colnames(mean_patch_size) <- c("id", "mean_area")

landscapes <- merge(landscapes, mean_patch_size, by = "id")

NP_PS_plot <- 
  ggplot(data = landscapes,
         mapping = aes(x = mean_area, y = number_patches)) +
  labs(x = "Mean patch size",
       y = "Number of forest patches") +
  geom_point() +
  geom_smooth(method = "lm", col = "black")+
  theme_pubr(base_size = 20)
NP_PS_plot

mod1 = lm(landscapes$mean_area~landscapes$number_patches)
summary(mod1)
cor.test(landscapes$mean_area,landscapes$number_patches)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")
ggsave("Mean_patch_size_vs_Number_of_patches.png", plot=NP_PS_plot,
       width = 200, height = 200, units = "mm",
       dpi = 600)

# Patch size evenness
calculate_pielou <- function(area_vector) {
  N <- sum(area_vector)  
  S <- length(area_vector)  
  p_i <- area_vector / N  
  shannon <- -sum(p_i * log(p_i))  
  pielou <- shannon / log(S)  
  return(pielou)
}

pielou_even <- aggregate(area ~ landscape_id, data = patch_area, FUN = calculate_pielou)

colnames(pielou_even) <- c("id", "evenness")
landscapes <- merge(landscapes, pielou_even, by = "id")

NP_PE_plot <- 
  ggplot(data = landscapes,
         mapping = aes(x = evenness, y = number_patches)) +
  labs(x = "Patch size evenness",
       y = "Number of forest patches") +
  geom_point() +
  geom_smooth(method = "lm", col = "black")+
  geom_hline(yintercept = c(20, 40))+
  theme_pubr(base_size = 20)
NP_PE_plot

cor.test(landscapes$evenness,landscapes$number_patches)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")
ggsave("Patch_size_evenness_vs_Number_of_patches.png", plot=NP_PE_plot,
       width = 200, height = 200, units = "mm",
       dpi = 600)


hist(log2(patch_area$area))
abline(v = c(-1, 2, 6), col = "red")

# Small 0.5 ha to 10 ha
# Medium > 10 ha to 50 ha
# Large > 50 ha
