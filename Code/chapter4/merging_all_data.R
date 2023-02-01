abundances <- read.csv("data/abundances.csv", sep = ";")
precip <- read.csv("data/precipitation_per_year.csv", sep = ",")
traits <- read.csv("data/plant_species_traits.csv", sep = ",")

nfd <- read.csv("data/NFD.csv", sep = ",")


library(dplyr)
library(ggplot2)
library(ellipse)
library(latex2exp)
library(ggpubr)


# Group the data by year, species, and compute the mean and sd of ND and FD
nd_fd <- nfd %>% 
  group_by(year, species) %>%
  summarize(ND = mean(ND, na.rm = TRUE),
            FD = mean(FD, na.rm = TRUE))
nd_fd$year <- as.factor(nd_fd$year)



ab <- abundances %>% 
  group_by(year, species) %>%
  summarize(indiv = sum(individuals, na.rm = TRUE))
ab$year <- as.factor(ab$year)



indiv <- NULL
weather <- NULL
func <- NULL
for (i in 1:nrow(nd_fd)){
  indiv <- c(indiv, ab$indiv[which(ab[which(ab$year == nd_fd$year[i]),]$species == nd_fd$species[i])])
  weather <- rbind(weather, precip[which(precip$year == nd_fd$year[i]),][, 2:11])
  func <- rbind(func, traits[which(traits$code == nd_fd$species[i]),][3:16])
}


nd_fd$indiv <- indiv
nd_fd <- cbind(nd_fd, weather)
nd_fd <- cbind(nd_fd, func)


write.csv(nd_fd, file = "Results/all_elements.csv", row.names = FALSE)




