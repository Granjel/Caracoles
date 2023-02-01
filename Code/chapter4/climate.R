# R squared
r2 <- function(glmodel){
  r <- 1 - (glmodel$deviance / glmodel$null.deviance)
  return(r)
}



data <- read.csv("Results/all_elements.csv")

library(dplyr)
library(ggplot2)
library(ellipse)
library(latex2exp)
library(ggpubr)
library(ggpmisc)




data <- data[complete.cases(data),]


d <- data[which(data$species %in% c("BETA", "BEMA", "CETE", "HOME", "LEMA", "PAIN", "POMA", "SASO")),]





nfd <- read.csv("data/NFD.csv", sep = ",")











summary(glm_ND_traits <- glm(ND ~ germination_rate + seed_survival + H + CS + LA + SLA + LAI + RD + SRL + TDMr + SRA + CN + C13 + N15, data = data))
summary(glm_FD_traits <- glm(FD ~ germination_rate + seed_survival + H + CS + LA + SLA + LAI + RD + SRL + TDMr + SRA + CN + C13 + N15, data = data))

summary(glm_ND_traits <- glm(ND ~ prec + variation + fall + spring + summer + winter + fall_rel + spring_rel + summer_rel + winter_rel, data = data))
summary(glm_FD_traits <- glm(FD ~ prec + fall + spring + summer + winter + fall_rel + spring_rel + summer_rel + winter_rel, data = data))

head(data)



library(vegan)
perma_ND_traits <- adonis2(data$ND ~ H + CS + LA + SLA + LAI + RD + SRL + TDMr + SRA + CN + C13 + N15, data = data, method = "euclidean")
perma_FD_traits <- adonis2(data$FD ~ H + CS + LA + SLA + LAI + RD + SRL + TDMr + SRA + CN + C13 + N15, data = data, method = "euclidean")

perma_ND_weather <- adonis2(data$ND ~ prec + fall + spring + summer + winter, data = data, method = "euclidean")
perma_FD_weather <- adonis2(data$FD ~ prec + fall + spring + summer + winter, data = data, method = "euclidean")


summary(glm(ND ~ winter, data = data))
summary(glm(ND ~ winter + winter^2, data = data))




#plants table
total_indiv <- data %>% group_by(species) %>% summarize(tot_indiv = sum(indiv))
total_indiv$tot_indiv <- (total_indiv$tot_indiv / sum(total_indiv$tot_indiv)) * 100

plants <- data.frame("CODE" = names(table(data$species)), "Species" = NA, "Years_number" = as.numeric(table(data$species)), "Individuals" = total_indiv$tot_indiv)

xtable::xtable(plants)




# Species' niches
# prec
ggplot(data = data, aes(x = prec, y = ND)) +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = "grey50", alpha = 0.15) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

data$prec2 <- data$prec^2
summary(glm(ND ~ prec, data = data))
summary(glm(ND ~ prec + prec2, data = data))



#87CEFA winter

seasons <- c("dodgerblue3", "forestgreen", "gold3", "tan4")

# winter
w_ND <-
  ggplot(data = data, aes(x = winter, y = ND)) +
  geom_point(color = seasons[1]) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = seasons[1], fill = seasons[1], alpha = 0.15) +
  scale_y_continuous(limits = c(1, 4), breaks = seq(1, 4, 0.5)) +
  labs(x = "Winter precipitation (mm)", y = "Species' niche") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

#data$winter2 <- data$winter^2
#summary(glm(ND ~ winter, data = data))
#summary(glm(ND ~ winter + winter2, data = data))
#r2(glm(ND ~ winter + winter2, data = data))


# fall
f_ND <-
  ggplot(data = data, aes(x = fall, y = ND)) +
  geom_point(color = seasons[4]) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = seasons[4], fill = seasons[4], alpha = 0.15) +
  scale_y_continuous(limits = c(1, 4), breaks = seq(1, 4, 0.5)) +
  labs(x = "Fall precipitation (mm)", y = "Species' niche") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.title.y = element_text(color = "white"))

#data$fall2 <- data$fall^2
#summary(glm(ND ~ fall, data = data))
#summary(glm(ND ~ fall + fall2, data = data))
#r2(glm(ND ~ fall + fall2, data = data))


# summer
s_ND <-
  ggplot(data = data, aes(x = summer, y = ND)) +
  geom_point(color = seasons[3]) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = seasons[3], fill = seasons[3], alpha = 0.15) +
  labs(x = "Summer precipitation (mm)", y = "Species' niche") +
  scale_y_continuous(limits = c(1, 4), breaks = seq(1, 4, 0.5)) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

#data$summer2 <- data$summer^2
#summary(glm(ND ~ summer, data = data))
#summary(glm(ND ~ summer + summer2, data = data))
#r2(glm(ND ~ summer + summer2, data = data))


# spring
g_ND <-
  ggplot(data = data, aes(x = spring, y = ND)) +
  geom_point(color = seasons[2]) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = seasons[2], fill = seasons[2], alpha = 0.15) +
  labs(x = "Spring precipitation (mm)", y = "Species' niche") +
  scale_y_continuous(limits = c(1, 4), breaks = seq(1, 4, 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.title.y = element_text(color = "white"))

#data$spring2 <- data$spring^2
#summary(glm(ND ~ spring, data = data))
#spring_ND <- glm(ND ~ spring + spring2, data = data)
#r2(spring_ND)


ggarrange(w_ND, g_ND, s_ND, f_ND, labels = "auto", align = "hv", font.label = list(size = 20))


ggsave(last_plot(), filename = "Figures/chapter4/rain_ND_seasons.png", device = "png", dpi = 320,
       width = 8, height = 8, units = "in")










# prec
t_FD <-
  ggplot(data = data, aes(x = prec, y = FD)) +
  geom_point(color = "black") +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = "grey50", fill = "grey50", alpha = 0.15) +
  labs(x = "Total precipitation (mm)", y = "Species' fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

ggsave(last_plot(), filename = "Figures/chapter4/rain_FD.png", device = "png", dpi = 320,
       width = 7, height = 6, units = "in")

data$prec2 <- data$prec^2
r2(glm(FD ~ prec, data = data))
prec_FD <- glm(FD ~ prec + prec2, data = data)
summary(prec_FD)
r2(prec_FD)


# winter
ggplot(data = data, aes(x = winter, y = FD)) +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = "grey50", alpha = 0.15) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

data$winter2 <- data$winter^2
summary(glm(FD ~ winter, data = data))
winter_FD <- glm(FD ~ winter + winter2, data = data)
summary(winter_FD)
r2(winter_FD)


# fall
ggplot(data = data, aes(x = fall, y = FD)) +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = "grey50", alpha = 0.15) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

data$fall2 <- data$fall^2
summary(glm(FD ~ fall, data = data))
fall_FD <- glm(FD ~ fall + fall2, data = data)
summary(fall_FD)
r2(fall_FD)


# summer
ggplot(data = data, aes(x = summer, y = FD)) +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = "grey50", alpha = 0.15) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

data$summer2 <- data$summer^2
summary(glm(FD ~ summer, data = data))
summer_FD <- glm(FD ~ summer + summer2, data = data)
summary(summer_FD)
r2(summer_FD)


# spring
ggplot(data = data, aes(x = spring, y = FD)) +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = "grey50", alpha = 0.15) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

data$spring2 <- data$spring^2
summary(glm(FD ~ spring, data = data))
spring_FD <- glm(FD ~ spring + spring2, data = data)
summary(spring_FD)
r2(spring_FD)






seasons <- c("dodgerblue3", "forestgreen", "gold3", "tan4")

# winter
w_FD <-
  ggplot(data = data, aes(x = winter, y = FD)) +
  geom_point(color = seasons[1]) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = seasons[1], fill = seasons[1], alpha = 0.15) +
  labs(x = "Winter precipitation (mm)", y = "Species' fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

#data$winter2 <- data$winter^2
#summary(glm(ND ~ winter, data = data))
#summary(glm(ND ~ winter + winter2, data = data))
#r2(glm(ND ~ winter + winter2, data = data))


# fall
f_FD <-
  ggplot(data = data, aes(x = fall, y = FD)) +
  geom_point(color = seasons[4]) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = seasons[4], fill = seasons[4], alpha = 0.15) +
  labs(x = "Fall precipitation (mm)", y = "Species' fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.title.y = element_text(color = "white"))

#data$fall2 <- data$fall^2
#summary(glm(ND ~ fall, data = data))
#summary(glm(ND ~ fall + fall2, data = data))
#r2(glm(ND ~ fall + fall2, data = data))


# summer
s_FD <-
  ggplot(data = data, aes(x = summer, y = FD)) +
  geom_point(color = seasons[3]) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = seasons[3], fill = seasons[3], alpha = 0.15) +
  labs(x = "Summer precipitation (mm)", y = "Species' fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank())

#data$summer2 <- data$summer^2
#summary(glm(ND ~ summer, data = data))
#summary(glm(ND ~ summer + summer2, data = data))
#r2(glm(ND ~ summer + summer2, data = data))


# spring
g_FD <-
  ggplot(data = data, aes(x = spring, y = FD)) +
  geom_point(color = seasons[2]) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = seasons[2], fill = seasons[2], alpha = 0.15) +
  labs(x = "Spring precipitation (mm)", y = "Species' fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.title.y = element_text(color = "white"))

#data$spring2 <- data$spring^2
#summary(glm(ND ~ spring, data = data))
#spring_ND <- glm(ND ~ spring + spring2, data = data)
#r2(spring_ND)


ggarrange(w_FD, g_FD, s_FD, f_FD, labels = "auto", align = "hv", font.label = list(size = 20))


ggsave(last_plot(), filename = "Figures/chapter4/rain_FD_seasons.png", device = "png", dpi = 320,
       width = 8, height = 8, units = "in")

























































