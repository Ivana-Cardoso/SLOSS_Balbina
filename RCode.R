# SLOSS cube hypothesis test
# in Balbina Hidroeletric Dam landscapes
#
# Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: February 21, 2024

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
new_crs <- CRS("+proj=utm +south +zone=21 +datum=WGS84 +units=m +no_defs")
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
proj4string(extent_grid) <- crs(new_crs)
points <- sp::makegrid(extent_grid, cellsize = 4000) # evenly spaced points (4km)
points(points, pch=".")


# Calculate number of forest patches in 2 km-radius landscape
# https://r-spatialecology.github.io/landscapemetrics/articles/guide_sample_lsm.html
grid <- st_as_sf(points, coords = c("x1", "x2"), crs = new_crs)
number_patches = landscapemetrics::sample_lsm(Forest_formation, y = grid, size=2000,
                                              shape = "square",
                                              what = "lsm_c_np",
                                              directions = 8)
number_patches = subset(number_patches, number_patches$class == 1)


# Calculate forest cover in 2 km-radius landscape
forest_cover = landscapemetrics::sample_lsm(Forest_formation, y = grid, size=2000,
                                            shape = "square",
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
  geom_hline(yintercept = c(20, 50))+
  theme_pubr(base_size = 20)
  
ggsave("Number_patches_vs_Forest_cover.png", plot=NP_FC_plot,
       width = 200, height = 200, units = "mm",
       dpi = 600)


# Select landscapes with 20 to 40% forest cover because, within these values, forest cover is visually not correlated with number of forest fragments
metrics = subset(metrics, metrics$forest_cover <= 50 & metrics$forest_cover >= 20)
cor(metrics$number_patches, metrics$forest_cover) # not correlated (r=-0.32)
selected_landscapes <- subset(metrics, 
                              metrics$id == 145 | metrics$id == 583 | 
                              metrics$id == 503 | metrics$id == 106 |
                              metrics$id == 112 | metrics$id == 428 |
                              metrics$id == 259 | metrics$id == 468 |
                              metrics$id == 143 | metrics$id == 229 |
                              metrics$id == 422 | metrics$id == 381 |
                              metrics$id == 188 | metrics$id == 228 |
                              metrics$id == 105 | metrics$id == 28 |
                              metrics$id == 226 | metrics$id == 186 |
                              metrics$id == 376 | metrics$id == 128 |
                              metrics$id == 624 | metrics$id == 545 |
                              metrics$id == 585 | metrics$id == 586 |
                              metrics$id == 546 | metrics$id == 76 |
                              metrics$id == 299 | metrics$id == 588 |
                              metrics$id == 387 | metrics$id == 427)

descarted_landscapes <- subset(metrics, metrics$id == 690 | metrics$id == 727 | 
                                metrics$id == 646 | metrics$id == 654 |
                                 metrics$id == 789 | metrics$id == 728 |
                                 metrics$id == 489 | metrics$id == 331 |
                                 metrics$id == 749 | metrics$id == 693 |
                                 metrics$id == 687 | metrics$id == 614 |
                                 metrics$id == 452)

NP_FC_plot +
geom_point(selected_landscapes, mapping = aes(x=number_patches, y=forest_cover), col = "red")+
geom_point(descarted_landscapes, mapping = aes(x=number_patches, y=forest_cover), col = "blue")


paired_landscapes <- read.csv("paired_landscapes.csv")
paired_landscapes$group <- as.factor(paired_landscapes$group)

paired = 
  ggplot(data = paired_landscapes,
       mapping = aes(x = number_patches, y = forest_cover,
                     group = group)) +
  labs(x = "Number of forest patches",
       y = "Forest cover (%)") +
  geom_point(size = 3, shape = 16) +
  geom_line()+
  theme_pubr(base_size = 20)

paired +
  geom_point(data = paired_landscapes[21:24,],
             mapping = aes(x = number_patches, y = forest_cover),
             color = "red")

correl = paired_landscapes[1:20,]

cor(correl$number_patches, correl$forest_cover)



# Export data to visualize in QGIS
plot(Forest_formation)
points(points, pch=".")
points$id = seq(1,nrow(points), 1)
landscapes = st_buffer(grid, endCapStyle = "SQUARE", dist = 2000)
landscapes$id = points$id
landscapes = landscapes[landscapes$id %in% metrics$id,]
landscapes = merge(landscapes, metrics, by="id")
plot(landscapes, add=T)

writeRaster(Forest_formation, "Floresta_Balbina", format="GTiff")
st_write(landscapes, "landscapes.shp", append = F)
write.csv(selected_landscapes, "selected_landscapes.csv")
