# Ellipses per year!

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


nd_fd$core <- "No"
nd_fd$core[which(nd_fd$species %in% names(which(table(nd_fd$species) == 5)))] <- "Yes"
nd_fd$core <- as.factor(nd_fd$core)



# Plot ellipses around each year's datapoints in a ND-FD scatterplot

ggplot(nd_fd, aes(x = ND, y = FD, color = factor(year))) + 
  geom_point(aes(shape = core), size = 3) +
  stat_ellipse(geom = "polygon", aes(fill = factor(year)), level = 0.95, alpha = 0.3, type = "t") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Species' niche", y = "Species' fitness", fill = "Year", color = "Year") +
  scale_fill_manual(values = c("#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7")) +
  scale_color_manual(values = c("#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7")) +
  scale_shape_manual(name = "Core\nspecies?", values = c(1, 16)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"))

ggsave(last_plot(), filename = "Figures/chapter4/ND_FD_year_core.png", device = "png", dpi = 320,
       width = 7, height = 6, units = "in")


ggplot(nd_fd, aes(x = ND, y = FD, color = factor(year))) + 
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", aes(fill = factor(year)), level = 0.95,alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Niche differences", y = "Fitness differences", fill = "Year", color = "Year") +
  scale_fill_manual(values = c("#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7")) +
  scale_color_manual(values = c("#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"))


ggsave(last_plot(), filename = "Figures/chapter4/ND_FD_year.png", device = "png", dpi = 320,
       width = 7, height = 6, units = "in")



#Pairwise t-test!
summary(glm(ND ~ year, data = nd_fd))
summary(glm(FD ~ year, data = nd_fd))

res.snd <- glm(FD ~ year, data = nd_fd)
pairwise.t.test(nd_fd$ND, nd_fd$year,
                p.adjust.method = "BH")
pairwise.t.test(nd_fd$FD, nd_fd$year,
                p.adjust.method = "BH")

















