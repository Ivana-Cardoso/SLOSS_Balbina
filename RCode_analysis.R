# SLOSS cube hypothesis test
# in Balbina Hydroelectric Dam landscapes
#
# Author: Ivana Cardoso - ivanawaters@gmail.com
#
# Created: September 09, 2025
# Last modification: September 10, 2025

library(dplyr)
library(car)
library(ggplot2)

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

comm_landscape <- complete %>%
  group_by(landscape) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

comm_landscape <- comm_landscape[,-c(2:9)]

comm_landscape <- 
  comm_landscape %>% 
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))

landscapes$richness <- rowSums(comm_landscape[,2:171]) # richness by landscape

landscapes$type <- NA
landscapes$type[1:8] <- "FL"
landscapes$type[9:12] <- "MS"
landscapes$type[13:20] <- "SS"
landscapes$type <- as.factor(landscapes$type)



# Are there differences in ant species richness (FL vs MS vs SS)? 
# ANOVA
mod <- aov(richness~type, data = landscapes)
summary(mod) # F-value = 1.95, p-value = 0.17
TukeyHSD(mod) # There are NO differences

shapiro.test(mod$residuals) # normal residuals (p=0.81)
leveneTest(mod) # variances are equal (p=0.64)

# Plotting the graph
plot1 <-
ggplot(data = landscapes,
       mapping = aes(x = type, y = richness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2.5, alpha = 0.6) +
  xlab("Landscape type") +
  ylab("Number of ant species") +
  
  theme_bw(base_size = 10) +
  
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25)) +
  
  annotate("text", x = 3, y = 53,
           hjust = 0, vjust = 0, size = 2.5,
           parse = T, label = as.character(expression(italic(F)*"-value = 1.95"))) +
  
  annotate("text", x = 3, y = 53,
           hjust = 0, vjust = 2, size = 2.5,
           parse = T, label = as.character(expression(italic(p)*"-value = 0.17")))

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/figures")
ggsave(
  plot = plot1,
  filename = "Richness_LandscapeType.png", dpi = 500,
  width = 8 * 1.5, height = 8, units = 'cm')   
  


# Are there differences in ant species richness (Richness vs Number of patches)? 
# Linear regression
mod2 <- lm(richness~number_patches, data = landscapes) 
summary(mod2) # adj-r2 = 0.20, slope = -0.21, p = 0.03

shapiro.test(mod2$residuals) # normal residuals (p=0.90)

# Plotting the graph
plot2 <-
  ggplot(data = landscapes,
         mapping = aes(x = number_patches, y = richness)) +
  
  geom_smooth(method = "lm", color = "black") +
  geom_point(size = 2.5, alpha = 0.6) +
  
  labs(x = "Number of patches in the landscape",
       y = "Number of ant species") +
  
  theme_bw(base_size = 10) +
  
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25)) +
  
  annotate("text", x = 50, y = 51,
           hjust = 0, vjust = 0, size = 2.5,
           parse = T,
           label = as.character(expression(italic(r)^{2}*""[adj]*" = 0.20"))) +
  
  annotate("text", x = 50, y = 51,
           hjust = 0, vjust = 2, size = 2.5,
           parse = T,
           label = as.character(expression(italic()*"slope = -0.21"))) +
  
  annotate("text", x = 50, y = 51,
           hjust = 0, vjust = 4, size = 2.5,
           parse = T, label = as.character(expression(italic(p)*"-value = 0.03")))

plot2

setwd("C:/Users/ivana/OneDrive/PhD_INPA/2.SLOSS_question/Analises/SLOSS_Balbina/figures")
ggsave(
  plot = plot2,
  filename = "Richness_NumberPatches.png", dpi = 500,
  width = 8 * 1.5, height = 8, units = 'cm')   
