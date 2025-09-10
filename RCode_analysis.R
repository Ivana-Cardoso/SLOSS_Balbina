# SLOSS cube hypothesis test
# in Balbina Hydroelectric Dam landscapes
#
# Author: Ivana Cardoso - ivanawaters@gmail.com
#
# Created: September 09, 2025
# Last modification: September 09, 2025

library(dplyr)

# Importing data
rm(list = ls()); gc()
set.seed(13)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/data")
comm <- read.csv("ants_comm.csv", header = T)
env <- read.csv("ants_env.csv", header = T)

complete <- merge(comm, env, by = "point_ID")
complete <- complete[,c(1,171:180,2:170)]

# Removing 100% FC landscapes
complete <- complete[-c(1:40),]

# Testing the correlation between forest cover and number of patches 
# of the sampled landscapes
landscapes <- complete[!duplicated(complete$landscape),]
landscapes <- landscapes[-c(1,5:180)]
cor.test(landscapes$number_patches, landscapes$forest_cover) 
# cor = -0.22, not correlated

# Gamma diversity - Total species richness
rownames(complete) <- complete$point_ID
complete$richness <- rowSums(complete[,c(12:180)]) # richness by point
complete$landscape <- as.factor(complete$landscape)

colSums(matrix(complete$richness, nrow = 10)) # essa nao é a solucao pq ta somando especies que ocorreram em pontos diferentes na paisagem. Deve contar so como 1 para a paisagem, então vou ter que corrigir.


