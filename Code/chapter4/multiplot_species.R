#Multiplot species

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



leg <- 
  ggplot(nd_fd, aes(x = year, y = FD)) + 
  geom_point(aes(color = species, shape = species), size = 3) + 
  scale_shape_manual(values = c(16, 16, 21, 16, 16, 16, 21, 21, 21, 16, 21, 21, 21, 16, 16, 21, 21, 21)) +
  labs(x = "Year", y = "Species' fitness", shape = "Species", color = "Species") +
  guides(color = guide_legend(ncol = 3), shape = guide_legend(ncol = 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "white"),
        axis.ticks = element_line(color = "white"),
        axis.line = element_line(color = "white"),
        axis.title = element_text(color = "white"),
        axis.text = element_text(color = "white"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black")) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "white") +
  theme(legend.position = c(0.4, 0.4))



p1 <- 
  ggplot(nd_fd, aes(x = ND, y = FD)) + 
  geom_point(aes(color = factor(species), shape = factor(core)), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Species' niche", y = "Species' fitness", shape = "Species", color = "Species") +
  scale_shape_manual(name = "Core\nspecies?", values = c(1, 16)) +
  scale_x_continuous(limits = c(1, 4), breaks = seq(1, 4, 0.5)) +
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "white"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")

#ggsave(p1, filename = "Figures/chapter4/ND_FD_species.png", device = "png", dpi = 320,
#       width = 8, height = 6, units = "in")


# FD and YEAR, species
p2 <- 
  ggplot(nd_fd, aes(x = year, y = FD)) + 
  geom_point(aes(color = factor(species), shape = factor(core)), size = 3) +
  #geom_point(aes(color = factor(species), shape = factor(species)), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(data = subset(nd_fd, core == "Yes"), aes(x = year, y = FD, group = species, color = species), alpha = 0.3) +
  labs(x = "Year", y = "Species' fitness", shape = "Species", color = "Species") +
  scale_shape_manual(name = "Core\nspecies?", values = c(1, 16)) +
  #scale_shape_manual(values = rep(15:18, len = 18)) +
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "white"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")

#ggsave(p2, filename = "Figures/chapter4/FD_year.png", device = "png", dpi = 320,
#       width = 8, height = 6, units = "in")


# ND and YEAR, species

#  ggplot(nd_fd, aes(y = year, x = ND)) + 
#  geom_point(aes(color = factor(species), shape = factor(species)), size = 3) +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#  labs(y = "Year", x = "Species' niche", shape = "Species", color = "Species") +
#  scale_shape_manual(values = rep(15:18, len = 18)) +
#  scale_x_continuous(limits = c(1, 4), breaks = seq(1, 4, 0.5)) +
#  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
#  theme_bw() +
#  theme(panel.grid = element_blank(),
#        text = element_text(size = 15, color = "black"),
#        axis.ticks = element_line(color = "black"),
#        legend.title.align = 0.5, legend.position = "none")

p3 <- 
  ggplot(nd_fd, aes(x = year, y = ND)) + 
  geom_point(aes(color = factor(species), shape = factor(core)), size = 3) +
  geom_line(data = subset(nd_fd, core == "Yes"), aes(x = year, y = ND, group = species, color = species), alpha = 0.3) +
  labs(x = "Year", y = "Species' niche", color = "Species") +
  scale_y_continuous(limits = c(1, 4), breaks = seq(1, 4, 0.5)) +
  scale_shape_manual(name = "Core\nspecies?", values = c(1, 16)) +
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  coord_flip()

#ggsave(p3, filename = "Figures/chapter4/ND_year.png", device = "png", dpi = 320,
#       width = 8, height = 6, units = "in")


ggarrange(p1, p2, p3, leg, labels = c("a", "b", "c", ""), common.legend = FALSE, align = "hv", font.label = list(size = 20))

ggsave(last_plot(), filename = "Figures/chapter4/species_coren.png", device = "png", dpi = 320,
       width = 8.5, height = 8, units = "in")












