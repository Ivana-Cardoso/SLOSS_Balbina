# SLOSS cube hypothesis test
# in Balbina Hydroelectric Dam landscapes
#
# Author: Ivana Cardoso - ivanawaters@gmail.com
#
# Last modification: March 28, 2024
#
# Select sample points in landscapes with LOW EVENNESS
#
# In landscapes with low evenness, I selected 4 points within the patches representing the landscape. For instance, if the landscape comprises small patches, then 4 of these patches were chosen. Additionally, I selected 3 patches from each of the other categories of patches

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
# Landscape 259
L259 <- landscapes[landscapes$id == 259,]
forest_landscapes <- crop(forest_formation, L259)
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
L259 <- extract_lsm(forest_landscapes, 
                   y = coords, 
                   extract_id = coords$id,
                   what = "lsm_p_area", 
                   directions = 8)
L259 <- L259[, c(6,7)]
colnames(L259) <- c("area", "id")
coords <- merge(coords, L259, by = "id")
coords$log_area <- log2(coords$area)

# Plot log2(area) against sampling points
plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

# Select patches with area from 0.5 to 4 ha
L259_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L259_S$landscape_id <- rep(259, length(L259_S$id))
L259_S$landscape <- rep("SS", length(L259_S$id))

L259_M = coords[coords$area > 10 & coords$area <= 50,]
L259_M$landscape_id <- rep(259, length(L259_M$id))
L259_M$landscape <- rep("SS", length(L259_M$id))

L259_L = coords[coords$area > 50,]
L259_L$landscape_id <- rep(259, length(L259_L$id))
L259_L$landscape <- rep("SS", length(L259_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L259_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L259_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L259_L$geometry)
patch_L <- forest[points_L]

sampling_points <- data.frame()

# Extract the values from the extent of landscape
landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

# Create points equally spaced at 100 meters within the extent of landscape i  
grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

# Select only the points in the patches
grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

# Converts spatial points object into a data.frame and identifies each one.  
grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

# Select random points at least 200 meters apart  
set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                           dist.min = 200, coord.X = "X", 
                           coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

# Converts the results into a data.frame, identifies, and adds the coordinates
samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(259, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(259, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(259, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- samp_points
write_sf(L259_S, "L259_S.shp")
write_sf(L259_M, "L259_M.shp")
write_sf(L259_L, "L259_L.shp")

# Landscape 301
L301 <- landscapes[landscapes$id == 301,]
forest_landscapes <- crop(forest_formation, L301)
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

L301 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L301 <- L301[, c(6,7)]
colnames(L301) <- c("area", "id")
coords <- merge(coords, L301, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L301_S = coords[coords$area >= 0.5  & coords$area <= 10,]
L301_S$landscape_id <- rep(301, length(L301_S$id))
L301_S$landscape <- rep("SS", length(L301_S$id))

L301_M = coords[coords$area > 10 & coords$area <= 50,]
L301_M$landscape_id <- rep(301, length(L301_M$id))
L301_M$landscape <- rep("SS", length(L301_M$id))

L301_L = coords[coords$area > 50,]
L301_L$landscape_id <- rep(301, length(L301_L$id))
L301_L$landscape <- rep("SS", length(L301_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L301_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L301_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L301_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(301, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(301, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(301, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L301_S, "L301_S.shp")
write_sf(L301_M, "L301_M.shp")
write_sf(L301_L, "L301_L.shp")

# Landscape 338
L338 <- landscapes[landscapes$id == 338,]
forest_landscapes <- crop(forest_formation, L338)
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

L338 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L338 <- L338[, c(6,7)]
colnames(L338) <- c("area", "id")
coords <- merge(coords, L338, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5,10,50))

L338_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L338_S$landscape_id <- rep(338, length(L338_S$id))
L338_S$landscape <- rep("SS", length(L338_S$id))

L338_M = coords[coords$area > 10 & coords$area <= 50,]
L338_M$landscape_id <- rep(338, length(L338_M$id))
L338_M$landscape <- rep("SS", length(L338_M$id))

L338_L = coords[coords$area > 50,]
L338_L$landscape_id <- rep(338, length(L338_L$id))
L338_L$landscape <- rep("SS", length(L338_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L338_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L338_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L338_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(338, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(338, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(338, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L338_S, "L338_S.shp")
write_sf(L338_M, "L338_M.shp")
write_sf(L338_L, "L338_L .shp")

# Landscape 468
L468 <- landscapes[landscapes$id == 468,]
forest_landscapes <- crop(forest_formation, L468)
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

L468 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L468 <- L468[, c(6,7)]
colnames(L468) <- c("area", "id")
coords <- merge(coords, L468, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L468_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L468_S$landscape_id <- rep(468, length(L468_S$id))
L468_S$landscape <- rep("SS", length(L468_S$id))

L468_M = coords[coords$area > 10 & coords$area <= 50,]
L468_M$landscape_id <- rep(468, length(L468_M$id))
L468_M$landscape <- rep("SS", length(L468_M$id))

L468_L = coords[coords$area > 50,]
L468_L$landscape_id <- rep(468, length(L468_L$id))
L468_L$landscape <- rep("SS", length(L468_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L468_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L468_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L468_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(468, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(468, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(468, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L468_S, "L_468_S.shp")
write_sf(L468_M, "L_468_M.shp")
write_sf(L468_L, "L_468_L.shp")

# Landscape 588
L588 <- landscapes[landscapes$id == 588,]
forest_landscapes <- crop(forest_formation, L588)
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

L588 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L588 <- L588[, c(6,7)]
colnames(L588) <- c("area", "id")
coords <- merge(coords, L588, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L588_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L588_S$landscape_id <- rep(588, length(L588_S$id))
L588_S$landscape <- rep("SS", length(L588_S$id))

L588_M = coords[coords$area > 10 & coords$area <= 50,]
L588_M$landscape_id <- rep(588, length(L588_M$id))
L588_M$landscape <- rep("SS", length(L588_M$id))

L588_L = coords[coords$area > 50,]
L588_L$landscape_id <- rep(588, length(L588_L$id))
L588_L$landscape <- rep("SS", length(L588_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L588_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L588_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L588_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(588, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(588, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(588, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L588_S, "L588_S.shp")
write_sf(L588_M, "L588_M.shp")
write_sf(L588_L, "L588_L.shp")

## MEDIUM SIZE ##
# Landscape 66
L66 <- landscapes[landscapes$id == 66,]
forest_landscapes <- crop(forest_formation, L66)
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

L66 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L66 <- L66[, c(6,7)]
colnames(L66) <- c("area", "id")
coords <- merge(coords, L66, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L66_S = coords[coords$area >= 0.5  & coords$area <= 10,]
L66_S$landscape_id <- rep(66, length(L66_S$id))
L66_S$landscape <- rep("MS", length(L66_S$id))

L66_M = coords[coords$area > 10 & coords$area <= 50,]
L66_M$landscape_id <- rep(66, length(L66_M$id))
L66_M$landscape <- rep("MS", length(L66_M$id))

L66_L = coords[coords$area > 50,]
L66_L$landscape_id <- rep(66, length(L66_L$id))
L66_L$landscape <- rep("MS", length(L66_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L66_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L66_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L66_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(66, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(66, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(66, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L66_S, "L66_S.shp")
write_sf(L66_M, "L66_M.shp")
write_sf(L66_L, "L66_L.shp")

# Landscape 143
L143 <- landscapes[landscapes$id == 143,]
forest_landscapes <- crop(forest_formation, L143)
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

L143 <- extract_lsm(forest_landscapes, 
                   y = coords, 
                   extract_id = coords$id,
                   what = "lsm_p_area", 
                   directions = 8)
L143 <- L143[, c(6,7)]
colnames(L143) <- c("area", "id")
coords <- merge(coords, L143, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L143_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L143_S$landscape_id <- rep(143, length(L143_S$id))
L143_S$landscape <- rep("MS", length(L143_S$id))

L143_M = coords[coords$area > 10 & coords$area <= 50,]
L143_M$landscape_id <- rep(143, length(L143_M$id))
L143_M$landscape <- rep("MS", length(L143_M$id))

L143_L = coords[coords$area > 50,]
L143_L$landscape_id <- rep(143, length(L143_L$id))
L143_L$landscape <- rep("MS", length(L143_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L143_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L143_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L143_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(143, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(143, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(143, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L143_S, "L143_S.shp")
write_sf(L143_M, "L143_M.shp")
write_sf(L143_L, "L143_L.shp")

# Landscape 263
L263 <- landscapes[landscapes$id == 263,]
forest_landscapes <- crop(forest_formation, L263)
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

L263 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L263 <- L263[, c(6,7)]
colnames(L263) <- c("area", "id")
coords <- merge(coords, L263, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L263_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L263_S$landscape_id <- rep(263, length(L263_S$id))
L263_S$landscape <- rep("MS", length(L263_S$id))

L263_M = coords[coords$area > 10 & coords$area <= 50,]
L263_M$landscape_id <- rep(263, length(L263_M$id))
L263_M$landscape <- rep("MS", length(L263_M$id))

L263_L = coords[coords$area > 50,]
L263_L$landscape_id <- rep(263, length(L263_L$id))
L263_L$landscape <- rep("MS", length(L263_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L263_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L263_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L263_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(263, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(263, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(263, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L263_S, "L263_S.shp")
write_sf(L263_M, "L263_M.shp")
write_sf(L263_L, "L263_L.shp")

# Landscape 340
L340 <- landscapes[landscapes$id == 340,]
forest_landscapes <- crop(forest_formation, L340)
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

L340 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L340 <- L340[, c(6,7)]
colnames(L340) <- c("area", "id")
coords <- merge(coords, L340, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L340_S = coords[coords$area >= 0.5  & coords$area <= 10,]
L340_S$landscape_id <- rep(340, length(L340_S$id))
L340_S$landscape <- rep("MS", length(L340_S$id))

L340_M = coords[coords$area > 10 & coords$area <= 50,]
L340_M$landscape_id <- rep(340, length(L340_M$id))
L340_M$landscape <- rep("MS", length(L340_M$id))

L340_L = coords[coords$area > 50,]
L340_L$landscape_id <- rep(340, length(L340_L$id))
L340_L$landscape <- rep("MS", length(L340_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L340_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L340_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L340_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(340, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(340, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(340, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L340_S, "L340_S.shp")
write_sf(L340_M, "L340_M.shp")
write_sf(L340_L, "L340_L.shp")

# Landscape 506
L506 <- landscapes[landscapes$id == 506,]
forest_landscapes <- crop(forest_formation, L506)
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

L506 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L506 <- L506[, c(6,7)]
colnames(L506) <- c("area", "id")
coords <- merge(coords, L506, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L506_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L506_S$landscape_id <- rep(506, length(L506_S$id))
L506_S$landscape <- rep("MS", length(L506_S$id))

L506_M = coords[coords$area > 10 & coords$area <= 50,]
L506_M$landscape_id <- rep(506, length(L506_M$id))
L506_M$landscape <- rep("MS", length(L506_M$id))

L506_L = coords[coords$area > 50,]
L506_L$landscape_id <- rep(506, length(L506_L$id))
L506_L$landscape <- rep("MS", length(L506_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L506_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L506_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L506_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(506, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(506, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(506, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L506_S, "L506_S.shp")
write_sf(L506_M, "L506_M.shp")
write_sf(L506_L, "L506_L.shp")

## LARGE SIZE ##
# Landscape 381
L381 <- landscapes[landscapes$id == 381,]
forest_landscapes <- crop(forest_formation, L381)
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

L381 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L381 <- L381[, c(6,7)]
colnames(L381) <- c("area", "id")
coords <- merge(coords, L381, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L381_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L381_S$landscape_id <- rep(381, length(L381_S$id))
L381_S$landscape <- rep("MS", length(L381_S$id))

L381_M = coords[coords$area > 10 & coords$area <= 50,]
L381_M$landscape_id <- rep(381, length(L381_M$id))
L381_M$landscape <- rep("MS", length(L381_M$id))

L381_L = coords[coords$area > 50,]
L381_L$landscape_id <- rep(381, length(L381_L$id))
L381_L$landscape <- rep("MS", length(L381_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L381_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L381_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L381_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(381, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(381, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(381, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L381_S, "L381_S.shp")
write_sf(L381_M, "L381_M.shp")
write_sf(L381_L, "L381_L.shp")

# Landscape 384
L384 <- landscapes[landscapes$id == 384,]
forest_landscapes <- crop(forest_formation, L384)
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

L384 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L384 <- L384[, c(6,7)]
colnames(L384) <- c("area", "id")
coords <- merge(coords, L384, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L384_S = coords[coords$area >= 0.5  & coords$area <= 10,]
L384_S$landscape_id <- rep(384, length(L384_S$id))
L384_S$landscape <- rep("MS", length(L384_S$id))

L384_M = coords[coords$area > 10 & coords$area <= 50,]
L384_M$landscape_id <- rep(384, length(L384_M$id))
L384_M$landscape <- rep("MS", length(L384_M$id))

L384_L = coords[coords$area > 50,]
L384_L$landscape_id <- rep(384, length(L384_L$id))
L384_L$landscape <- rep("MS", length(L384_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L384_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L384_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L384_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(384, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(384, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(384, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L384_S, "L384_S.shp")
write_sf(L384_M, "L384_M.shp")
write_sf(L384_L, "L384_L.shp")

# Landscape 422
L422 <- landscapes[landscapes$id == 422,]
forest_landscapes <- crop(forest_formation, L422)
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

L422 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L422 <- L422[, c(6,7)]
colnames(L422) <- c("area", "id")
coords <- merge(coords, L422, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L422_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L422_S$landscape_id <- rep(422, length(L422_S$id))
L422_S$landscape <- rep("MS", length(L422_S$id))

L422_M = coords[coords$area > 10 & coords$area <= 50,]
L422_M$landscape_id <- rep(422, length(L422_M$id))
L422_M$landscape <- rep("MS", length(L422_M$id))

L422_L = coords[coords$area > 50,]
L422_L$landscape_id <- rep(422, length(L422_L$id))
L422_L$landscape <- rep("MS", length(L422_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L422_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L422_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L422_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(422, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(422, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(422, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L422_S, "L422_S.shp")
write_sf(L422_M, "L422_M.shp")
write_sf(L422_L, "L422_L.shp")

# Landscape 503
L503 <- landscapes[landscapes$id == 503,]
forest_landscapes <- crop(forest_formation, L503)
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

L503 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L503 <- L503[, c(6,7)]
colnames(L503) <- c("area", "id")
coords <- merge(coords, L503, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L503_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L503_S$landscape_id <- rep(503, length(L503_S$id))
L503_S$landscape <- rep("MS", length(L503_S$id))

L503_M = coords[coords$area > 10 & coords$area <= 50,]
L503_M$landscape_id <- rep(503, length(L503_M$id))
L503_M$landscape <- rep("MS", length(L503_M$id))

L503_L = coords[coords$area > 50,]
L503_L$landscape_id <- rep(503, length(L503_L$id))
L503_L$landscape <- rep("MS", length(L503_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L503_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L503_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L503_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_M <- select.site(x = grid_M, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_M <- as.data.frame(samp_points_M)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(503, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_M) <- "id"
samp_points_M <- merge(samp_points_M, grid_M, by = "id")
samp_points_M$landscape_id <- rep(503, length(samp_points_M$id))
samp_points_M$patch <- rep("M", length(samp_points_M$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(503, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_M, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L503_S, "L503_S.shp")
write_sf(L503_M, "L503_M.shp")
write_sf(L503_L, "L503_L.shp")

# Landscape 549
L549 <- landscapes[landscapes$id == 549,]
forest_landscapes <- crop(forest_formation, L549)
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

L549 <- extract_lsm(forest_landscapes, 
                    y = coords, 
                    extract_id = coords$id,
                    what = "lsm_p_area", 
                    directions = 8)
L549 <- L549[, c(6,7)]
colnames(L549) <- c("area", "id")
coords <- merge(coords, L549, by = "id")
coords$log_area <- log2(coords$area)

plot(sort(coords$area), 1:length(coords$area))
abline(v=c(0.5, 10, 50))

L549_S = coords[coords$area >= 0.5 & coords$area <= 10,]
L549_S$landscape_id <- rep(549, length(L549_S$id))
L549_S$landscape <- rep("MS", length(L549_S$id))

L549_M = coords[coords$area > 10 & coords$area <= 50,]
L549_M$landscape_id <- rep(549, length(L549_M$id))
L549_M$landscape <- rep("MS", length(L549_M$id))

L549_L = coords[coords$area > 50,]
L549_L$landscape_id <- rep(549, length(L549_L$id))
L549_L$landscape <- rep("MS", length(L549_L$id))

## Select sampling points ##
forest <- rasterToPolygons(forest_landscapes, fun=function(x){x==1}, dissolve = TRUE)
forest <- st_as_sfc(forest)
forest <- st_cast(forest, "POLYGON")
points_S <- st_as_sf(L549_S$geometry)
patch_S <- forest[points_S]
points_M <- st_as_sf(L549_M$geometry)
patch_M <- forest[points_M]
points_L <- st_as_sf(L549_L$geometry)
patch_L <- forest[points_L]

landscape <- as(extent(forest_landscapes), 'SpatialPolygons')
proj4string(landscape) <- crs(new_crs)

grid <- sp::makegrid(landscape, cellsize = 100)
grid <- st_as_sf(grid, coords = c("x1", "x2"), crs = new_crs)

grid_S <- grid[patch_S,]
grid_M <- grid[patch_M,]
grid_L <- grid[patch_L,]

grid_S <- st_coordinates(grid_S)
grid_S <- as.data.frame(grid_S)
grid_S$id <- 1:nrow(grid_S)

grid_M <- st_coordinates(grid_M)
grid_M <- as.data.frame(grid_M)
grid_M$id <- 1:nrow(grid_M)

grid_L <- st_coordinates(grid_L)
grid_L <- as.data.frame(grid_L)
grid_L$id <- 1:nrow(grid_L)

set.seed(13)
samp_points_S <- select.site(x = grid_S, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 4, Nsets = 1, Nmin = 4)
# There's no medium size patch
# samp_points_M <- select.site(x = grid_M, var.site = "id", 
#                             dist.min = 200, coord.X = "X", 
#                             coord.Y = "Y", Nmax = 3, Nsets = 1, Nmin = 3)
set.seed(13)
samp_points_L <- select.site(x = grid_L, var.site = "id", 
                             dist.min = 200, coord.X = "X", 
                             coord.Y = "Y", Nmax = 6, Nsets = 1, Nmin = 6)

samp_points_S <- as.data.frame(samp_points_S)
samp_points_L <- as.data.frame(samp_points_L)

colnames(samp_points_S) <- "id"
samp_points_S <- merge(samp_points_S, grid_S, by = "id")
samp_points_S$landscape_id <- rep(549, length(samp_points_S$id))
samp_points_S$patch <- rep("S", length(samp_points_S$id))

colnames(samp_points_L) <- "id"
samp_points_L <- merge(samp_points_L, grid_L, by = "id")
samp_points_L$landscape_id <- rep(549, length(samp_points_L$id))
samp_points_L$patch <- rep("L", length(samp_points_L$id))

samp_points <- rbind(samp_points_S, samp_points_L)
plot(forest_landscapes)
points(samp_points$X, samp_points$Y, pch = 16, col = "blue")
sampling_points <- rbind(sampling_points, samp_points)
write_sf(L549_S, "L549_S.shp")
write_sf(L549_M, "L549_M.shp")
write_sf(L549_L, "L549_L.shp")

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/Exported_files")
sampling_points <- st_as_sf(sampling_points, coords = c("X", "Y"))
write_sf(sampling_points, "points_ldscp_low_evenness.shp")

