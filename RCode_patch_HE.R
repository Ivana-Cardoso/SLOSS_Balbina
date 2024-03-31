# SLOSS cube hypothesis test
# in Balbina Hydroelectric Dam landscapes
#
# Author: Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: March 28, 2024
#
# Select sample points in landscapes with HIGH EVENNESS
#
# In landscapes with high evenness, I selected 10 points in the same patch size section

rm(list = ls()); gc(); dev.off()

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/Exported_files")

# Load packages
library(raster)
library(sf)
library(landscapemetrics)

# Import data
forest_formation <- raster("forest_formation.tif")
landscapes<- read_sf("landscapes.shp")
new_crs <- CRS("+proj=utm +south +zone=21 +datum=WGS84 +units=m +no_defs")

# Clear workspace
gc()

# I modified Pavel's function (with his permission) to randomly select only 10 points in each landscape while respecting a minimum distance of 200 meters. The original function is available at: https://github.com/pdodonov/sampling
setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina")
source("siteSelection_v0.8.R")


#### SELECT PATCHES ####

## SMALL SIZE ##
# Landscape 76
L76 <- landscapes[landscapes$id == 76,]
forest_landscapes <- crop(forest_formation, L76)
patches <- get_patches(forest_landscapes, class = 1)

# Get centroids of patches
coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

# Plot landscape and sampling points
plot(forest_landscapes)
points(coords$x, coords$y)

# Convert coordinates to spatial data frame
coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

# Extract landscape metrics for sampling points
L76 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L76 <- L76[, c(6,7)]
colnames(L76) <- c("area", "id")
coords <- merge(coords, L76, by = "id")
coords$log_area <- log2(coords$area)

# Plot area 
plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5,10,50))

# Select patches with area from 0.5 to 10 ha
L76 = coords[coords$area >=  0.5 & coords$area <= 10,]
L76$landscape_id <- rep(76, length(L76$id))
L76$landscape <- rep("SS", length(L76$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L76$geometry)
patch <- forest[points]
  
sampling_points <- data.frame()

# Extract the values from the extent of landscape
  landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
  proj4string(landscape) <- crs(new_crs)
  
# Create points equally spaced at 100 meters within the extent of landscape i  
  grid <- sp::makegrid(landscape, cellsize = 100)
  grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)
  
# Select only the points in the patches
  grid <- grid[patch,]
  
# Converts spatial points object into a data.frame and identifies each one.  
  grid <- st_coordinates(grid)
  grid <- as.data.frame(grid)
  grid$id <- 1:nrow(grid)
  
# Select random points at least 200 meters apart  
  set.seed(13)
  samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)
  
# Converts the results into a data.frame, identifies, and adds the coordinates
  samp_points <- as.data.frame(samp_points)
  colnames(samp_points) <- "id"
  samp_points <- merge(samp_points, grid, by = "id")
  samp_points$landscape_id <- rep(76, length(samp_points$id))
    
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- samp_points
write_sf(L76, "L76.shp") #Exporting all possible points

# Landscape 387
L387 <- landscapes[landscapes$id == 387, ]
forest_landscapes <- crop(forest_formation, L387)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L387 <- extract_lsm(forest_landscapes, 
                   y = coords, 
                   extract_id = coords$id,
                   what = "lsm_p_area", 
                   directions = 8)
L387 <- L387[, c(6,7)]
colnames(L387) <- c("area", "id")
coords <- merge(coords, L387, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5,10,50))

L387 = coords[coords$area >=  0.5 & coords$area <= 10,]
L387$landscape_id <- rep(387, length(L387$id))
L387$landscape <- rep("SS", length(L387$id))

df.patches <- rbind(L76, L387)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L387$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)
grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(387, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")

sampling_points <- rbind(sampling_points, samp_points)
write_sf(L387, "L387.shp") #Exporting all possible points

# Landscape 545
L545 <- landscapes[landscapes$id == 545, ]
forest_landscapes <- crop(forest_formation, L545)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L545 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L545 <- L545[, c(6,7)]
colnames(L545) <- c("area", "id")
coords <- merge(coords, L545, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5,10,50))

L545 = coords[coords$area >=  0.5 & coords$area <= 10,]
L545$landscape_id <- rep(545, length(L545$id))
L545$landscape <- rep("SS", length(L545$id))

df.patches <- rbind(df.patches, L545)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L545$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(545, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")

sampling_points <- rbind(sampling_points, samp_points)
write_sf(L545, "L545.shp") #Exporting all possible points

# Landscape 584
L584 <- landscapes[landscapes$id == 584, ]
forest_landscapes <- crop(forest_formation, L584)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L584 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L584 <- L584[, c(6,7)]
colnames(L584) <- c("area", "id")
coords <- merge(coords, L584, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5,10,50))

L584 = coords[coords$area >=  0.5 & coords$area <= 10,]
L584$landscape_id <- rep(584, length(L584$id))
L584$landscape <- rep("SS", length(L584$id))

df.patches <- rbind(df.patches, L584)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L584$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(584, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")

sampling_points <- rbind(sampling_points, samp_points)
write_sf(L584, "L584.shp") #Exporting all possible points

# Landscape 586
L586 <- landscapes[landscapes$id == 586, ]
forest_landscapes <- crop(forest_formation, L586)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L586 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L586 <- L586[, c(6,7)]
colnames(L586) <- c("area", "id")
coords <- merge(coords, L586, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5,10,50))

L586 = coords[coords$area >=  0.5 & coords$area <= 10,]
L586$landscape_id <- rep(586, length(L586$id))
L586$landscape <- rep("SS", length(L586$id))

df.patches <- rbind(df.patches, L586)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L586$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(586, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")

sampling_points <- rbind(sampling_points, samp_points)
write_sf(L586, "L586.shp")

## MEDIUM SIZE ##
# Landscape 74
L74 <- landscapes[landscapes$id == 74, ]
forest_landscapes <- crop(forest_formation, L74)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L74 <- extract_lsm(forest_landscapes, 
                   y = coords, 
                   extract_id = coords$id,
                   what = "lsm_p_area", 
                   directions = 8)
L74 <- L74[, c(6,7)]
colnames(L74) <- c("area", "id")
coords <- merge(coords, L74, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L74 = coords[coords$area > 10 & coords$area <= 50,]
L74$landscape_id <- rep(74, length(L74$id))
L74$landscape <- rep("MS", length(L74$id))

df.patches <- rbind(df.patches, L74)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L74$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(74, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")

sampling_points <- rbind(sampling_points, samp_points)
write_sf(L74, "L74.shp")

# Landscape 115
L115 <- landscapes[landscapes$id == 115, ]
forest_landscapes <- crop(forest_formation, L115)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L115 <- extract_lsm(forest_landscapes, 
                   y = coords, 
                   extract_id = coords$id,
                   what = "lsm_p_area", 
                   directions = 8)
L115 <- L115[, c(6,7)]
colnames(L115) <- c("area", "id")
coords <- merge(coords, L115, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L115 = coords[coords$area > 10 & coords$area <= 50,]
L115$landscape_id <- rep(115, length(L115$id))
L115$landscape <- rep("MS", length(L115$id))

df.patches <- rbind(df.patches, L115)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L115$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(115, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")

sampling_points <- rbind(sampling_points, samp_points)
write_sf(L115, "L115.shp")

# Landscape 189
L189 <- landscapes[landscapes$id == 189, ]
forest_landscapes <- crop(forest_formation, L189)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L189 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L189 <- L189[, c(6,7)]
colnames(L189) <- c("area", "id")
coords <- merge(coords, L189, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L189 = coords[coords$area > 10 & coords$area <= 50,]
L189$landscape_id <- rep(189, length(L189$id))
L189$landscape <- rep("MS", length(L189$id))

df.patches <- rbind(df.patches, L189)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L189$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(189, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")

sampling_points <- rbind(sampling_points, samp_points)
write_sf(L189, "L189.shp")

# Landscape 269
L269 <- landscapes[landscapes$id == 269, ]
forest_landscapes <- crop(forest_formation, L269)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L269 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L269 <- L269[, c(6,7)]
colnames(L269) <- c("area", "id")
coords <- merge(coords, L269, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L269 = coords[coords$area > 10 & coords$area <= 50,]
L269$landscape_id <- rep(269, length(L269$id))
L269$landscape <- rep("MS", length(L269$id))

df.patches <- rbind(df.patches, L269)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L269$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(269, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")

sampling_points <- rbind(sampling_points, samp_points)
write_sf(L269, "L269.shp")

# Landscape 307
L307 <- landscapes[landscapes$id == 307, ]
forest_landscapes <- crop(forest_formation, L307)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L307 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L307 <- L307[, c(6,7)]
colnames(L307) <- c("area", "id")
coords <- merge(coords, L307, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L307 = coords[coords$area > 10 & coords$area <= 50,]
L307$landscape_id <- rep(307, length(L307$id))
L307$landscape <- rep("MS", length(L307$id))

df.patches <- rbind(df.patches, L307)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L307$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(307, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L307, "L307.shp")

## LARGE SIZE ##
# Landscape 112
L112 <- landscapes[landscapes$id == 112, ]
forest_landscapes <- crop(forest_formation, L112)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L112 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L112 <- L112[, c(6,7)]
colnames(L112) <- c("area", "id")
coords <- merge(coords, L112, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L112 = coords[coords$area > 50,]
L112$landscape_id <- rep(112, length(L112$id))
L112$landscape <- rep("FL", length(L112$id))

df.patches <- rbind(df.patches, L112)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L112$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(112, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L112, "L112.shp")

# Landscape 145
L145 <- landscapes[landscapes$id == 145,]
forest_landscapes <- crop(forest_formation, L145)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L145 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L145 <- L145[, c(6,7)]
colnames(L145) <- c("area", "id")
coords <- merge(coords, L145, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L145 = coords[coords$area > 50,]
L145$landscape_id <- rep(145, length(L145$id))
L145$landscape <- rep("FL", length(L145$id))

df.patches <- rbind(df.patches, L145)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L145$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(145, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L145, "L145.shp")

# Landscape 150
L150 <- landscapes[landscapes$id == 150,]
forest_landscapes <- crop(forest_formation, L150)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L150 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L150 <- L150[, c(6,7)]
colnames(L150) <- c("area", "id")
coords <- merge(coords, L150, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L150 = coords[coords$area > 50,]
L150 = coords[coords$log_area > 6,]
L150$landscape_id <- rep(150, length(L150$id))
L150$landscape <- rep("FL", length(L150$id))

df.patches <- rbind(df.patches, L150)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L150$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(150, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L150, "L150.shp")

# Landscape 186
L186 <- landscapes[landscapes$id == 186,]
forest_landscapes <- crop(forest_formation, L186)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L186 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L186 <- L186[, c(6,7)]
colnames(L186) <- c("area", "id")
coords <- merge(coords, L186, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L186 = coords[coords$area > 50,]
L186$landscape_id <- rep(186, length(L186$id))
L186$landscape <- rep("FL", length(L186$id))

df.patches <- rbind(df.patches, L186)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L186$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(186, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L186, "L186.shp")

# Landscape 228
L228 <- landscapes[landscapes$id == 228,]
forest_landscapes <- crop(forest_formation, L228)
patches <- get_patches(forest_landscapes, class = 1)

coords <- get_centroids(patches$layer_1, 
                        directions = 8, 
                        cell_center = TRUE, 
                        return_vec = FALSE, 
                        verbose = TRUE)
coords <- coords[, c(4:6)]
coords$id <- 1:length(coords$id)

plot(forest_landscapes)
points(coords$x, coords$y)

coords <- st_as_sf(coords, 
                   coords = c("x", "y"), 
                   crs = new_crs)

L228 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L228 <- L228[, c(6,7)]
colnames(L228) <- c("area", "id")
coords <- merge(coords, L228, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L228 = coords[coords$area > 50,]
L228$landscape_id <- rep(228, length(L228$id))
L228$landscape <- rep("FL", length(L228$id))

df.patches <- rbind(df.patches, L228)

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points <- st_as_sf(L228$geometry)
patch <- forest[points]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid <- grid[patch,]

grid <- st_coordinates(grid)
grid <- as.data.frame(grid)
grid$id <- 1:nrow(grid)

set.seed(13)
samp_points <- select.site(x = grid, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 10, Nsets = 1, Nmin = 10)

samp_points <- as.data.frame(samp_points)
colnames(samp_points) <- "id"
samp_points <- merge(samp_points, grid, by = "id")
samp_points$landscape_id <- rep(228, length(samp_points$id))

plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L228, "L228.shp")

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/Exported_files")
sampling_points <- st_as_sf(sampling_points, coords = c("X", "Y"))
write_sf(sampling_points, "points_ldscp_high_evenness.shp")

