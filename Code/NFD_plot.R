#######################################################################################################----
# Package requirement 
#######################################################################################################----
# required packages
#install.packages("ggalt")
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}# to handle dataframe
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)} # to make cool graph
if(!require(grid)) {install.packages("grid"); library(grid)}
if(!require(ggpubr)) {install.packages("ggpubr"); library(ggpubr)}
if(!require(gtable)) {install.packages("gtable"); library(gtable)}
if(!require(psych)) {install.packages("psych"); library(psych)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(ggthemes)) {install.packages("ggthemes"); library(ggthemes)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(ggalt)) {install.packages("ggalt"); library(ggalt)} #  Extra Coordinate Systems, 'Geoms', Statistical Transformations, Scales and Fonts for 'ggplot2'
if(!require(plotly)) {install.packages("plotly"); library(plotly)} # to make intergratif graphs
if(!require(matrixStats)) {install.packages("matrixStats"); library(matrixStats)}
if(!require(data.table)) {install.packages("data.table"); library(data.table)}

#========================================================================================================
# Set Working directory
#========================================================================================================
setwd("/Users/lisabuche/Documents/Stage/Project/N-F-plane")


#========================================================================================================
#==========================
# Modify the data_frame
#==========================
#========================================================================================================

#######################################################################################################----
# Data_frame construction
#######################################################################################################----
# import data
NFD <- read.csv("data/NFD_uncert.csv")
#NFD <- read.csv("/Users/lisabuche/Documents/Stage/Project/data/NFD_data_nlminb.csv")

# make the table vertical instead or horyzontal
str(NFD)
names(NFD)
ND_factor <- names(NFD)[2:19]
FD_factor <- names(NFD)[20:37]
equi_factor <- names(NFD)[38:55]
levels(as.factor(NFD$case))


NFD <- gather(NFD,all_of(ND_factor), key='species.ND', value="ND")
NFD <- gather(NFD,all_of(FD_factor), key='species.FD', value="FD")
NFD <- gather(NFD,all_of(equi_factor), key='species.equi', value="equi")
NFD <- drop_na(NFD , c(ND,FD))
NFD <- separate(NFD, species.ND, into= c("toDelete","species"), sep = "_")
NFD <- separate(NFD, species.FD, into= c("toDelete","species"), sep = "_")
NFD <- separate(NFD, species.equi, into= c("toDelete","species"), sep = "_")
NFD <- separate(NFD, case, into= c("year","plot.number"), sep = "_")

NFD  <- subset(NFD, select=c("year","plot.number","ND","FD","species" ,"equi" ))

# change the FD with the modified definition
#NFD$FD <-  1-1/(1-NFD$FD)

# Remove very high values, otherwise they bend the mean to -inf or +inf


#Remove outliers
# Can diminish the restriction by changing probs
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

NFD$ND <- remove_outliers(NFD$ND)
NFD$FD <- remove_outliers(NFD$FD)
NFD <- remove_missing(NFD)
#check for outlier
ggplot(NFD,aes(ND,FD,group=species)) +geom_boxplot()
ggplot(NFD,aes(FD,ND,group=species)) +geom_boxplot()

write_csv(NFD,
          file.path("/Users/lisabuche/Documents/Stage/Project/N-F-plane/data",
                    "NFD.csv"),
          na = "NA", append = F,
          col_names =TRUE)

#######################################################################################################----
# Computation of centroids
#######################################################################################################----
#----compute centroide per class with class equal = species----
centroids <- aggregate(cbind(ND,FD)~species,NFD,median)
f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
#If you want to calculate, say, 95% confidence instead of std. error, replace
#f <- function(z) qt(0.025,df=length(z)-1, lower.tail=F)* sd(z)/sqrt(length(z)) 

se        <- aggregate(cbind(se.x=ND,se.y=FD)~species,NFD,f)
centroids_species <- merge(centroids,se,   by = "species") 

##---- or class = species_ year ----
NFD <- unite(NFD, species, year, col= "species_year", sep = "_", remove = F)
centroids <- aggregate(cbind(ND,FD) ~ species_year, NFD, median)
f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
#If you want to calculate, say, 95% confidence instead of std. error, replace
#f <- function(z) qt(0.025,df=length(z)-1, lower.tail=F)* sd(z)/sqrt(length(z)) 

se        <- aggregate(cbind(se.x=ND,se.y=FD)~species_year,NFD,f)
centroids.condensed <- merge(centroids,se, by="species_year")
centroids_sp_y <- separate(centroids.condensed, species_year, 
                           into=c("species","year"), sep="_", remove = F)

#######################################################################################################----
# Computation of Functional Dispersion - fdisp
#######################################################################################################----
nfd.dispersion.fun <- function(nd,fd){
  mean_nd <- rowMedians(nd,rows = NULL, cols = NULL, na.rm = T, dim.= c(1,length(nd)))
  mean_fd <- rowMedians(fd,rows = NULL, cols = NULL, na.rm = T, dim.= c(1,length(fd)))
  euclid_dist_nfd <- mean(sqrt(((mean_nd - nd)^2 + (mean_fd - fd)^2)))
  return(c(mean_nd,mean_fd, euclid_dist_nfd))
}

nfd.dispersion.df <- data.frame(matrix(ncol=5, nrow=0))
names(nfd.dispersion.df) <- c("species","year","nd","fd",
                              "dispersion")
nfd.dispersion.df$species <- as.character(nfd.dispersion.df$species)
nfd.dispersion.df$year <- as.character(nfd.dispersion.df$year)


nfd.dispersion.df_sp<- data.frame(matrix(ncol=4, nrow=0))
names(nfd.dispersion.df_sp) <- c("species","nd","fd",
                                 "dispersion")
nfd.dispersion.df_sp$species <- as.character(nfd.dispersion.df_sp$species)

for ( species.id in levels(as.factor(centroids_species$species))){   
  nfd.to.add_sp <- nfd.dispersion.fun(NFD$ND[ NFD$species==species.id ],
                                      NFD$FD[ NFD$species==species.id])
  nfd.dispersion.df_sp <- add_row(nfd.dispersion.df_sp, 
                                   species = species.id, 
                                   nd =   nfd.to.add_sp[1],
                                   fd =   nfd.to.add_sp[2],
                                   dispersion=   nfd.to.add_sp[3])
  
  for( year.id in levels(as.factor(centroids_sp_y$year))){
    nfd.to.add <- nfd.dispersion.fun(NFD$ND[ NFD$species==species.id & NFD$year==year.id ],
                                     NFD$FD[ NFD$species==species.id & NFD$year==year.id ])
    nfd.dispersion.df <- add_row(nfd.dispersion.df, 
                                 species=species.id,
                                 year= year.id, 
                                 nd = nfd.to.add[1],
                                 fd = nfd.to.add[2],
                                 dispersion= nfd.to.add[3]
                                 
    )
  }
}

str(nfd.dispersion.df)


#========================================================================================================
#==========================
# Graphiques
#==========================
#========================================================================================================


#######################################################################################################----
# plot of the centroids
#######################################################################################################----


#---- centroide per class with class equal = species----
NFD_point <- ggplot(NFD,aes(color=factor(species), shape= year,
                            ND,FD)) +
  geom_errorbar(data=centroids_sp_y ,aes(ymin=FD-se.y,ymax=FD+se.y),
                width=0.1, alpha= 0.6)+
  geom_errorbarh(data=centroids_sp_y,aes(xmin=ND-se.x,xmax=ND+se.x,height=.005),
                 alpha = 0.6) + 
  geom_point(data=centroids_sp_y ,size=2) + 
   #geom_encircle() +
   coord_cartesian(xlim = c(-1, 5), ylim = c(-0.1, 0.5)) +
  theme_bw() 
NFD_point
# to make it interactive
ggplotly(NFD_point)
ggsave( plot= NFD_point,
        filename = "NFD_point_sp.pdf",
        path = "Figures/")

# you ca save it as a web page to share NFD_point_species

#---- or class = species_ year ----

NFD_plane_circle <- ggplot(centroids_sp_y ,aes(ND,FD,color=species, shape= year)) +
  #geom_point(data=centroids,size=2) + 
  geom_encircle(aes(group= species, fill=species,alpha=0.4),
                s_shape=1, expand=0,#spread=0,
                stat = "identity", position = "identity") + 
  #coord_cartesian(xlim = c(-3, 12), ylim = c(-15, 0)) +
  theme_bw() 
NFD_plane_circle
ggsave( plot= NFD_plane_circle,
        filename = "NFD_plane_circle.pdf",
        path = "/Users/lisabuche/Documents/Stage/Project/N-F-plane/Figures")

centroids_sp_y_abund <- subset(centroids_sp_y, species== "CHFU" |
                                 species== "HOMA"|
                                 species== "LEMA"|
                                 species== "CETE"|
                                 species== "SPRU"
                               )
NFD_point_sp_y <- ggplot(centroids_sp_y_abund  ,
                         aes(ND,FD,color=species, shape= year)) +
  #geom_point(data=centroids,size=2) + 
  geom_point() +
  geom_encircle(aes(group= species, fill=species),alpha=0.26,
                s_shape=0.5, expand=0,#spread=0,
                stat = "identity", position = "identity") + 
  coord_cartesian(xlim = c(0, 5), ylim = c(-0.1, 0.4)) +
  theme_bw() 

NFD_point_sp_y
ggsave( plot= NFD_point_sp_y,
  filename = "NFD_point_sp_y.pdf",
  path = "Figures/")
# to make it interactive
ggplotly(NFD_point_sp_y)
# you ca save it as a web page to share NFD_point_species_yr

NFD_plane_year <- ggplot(centroids_sp_y ,aes(ND,FD,color=year)) +
  #geom_point(data=centroids,size=2) + 
  geom_point() + 
  geom_encircle(aes(group= year, fill=year),alpha=0.4,
                s_shape=0.5, expand=0,#spread=0, 
                # expand = to expand the circle beyond the points
                # s_shape = 1 -> triangle, <1 -> rounder
                stat = "identity", position = "identity") + 
  #coord_cartesian(xlim = c(-3, 12), ylim = c(-15, 0)) +
  theme_bw() 
NFD_plane_year
ggsave( plot= NFD_plane_year,
        filename = "NFD_plane_year.pdf",
        path = "Figures/")



#######################################################################################################----
# Functional Dispersion graphs
#######################################################################################################----

nfd.dispersion.plot_sp <- ggplot(nfd.dispersion.df_sp , 
                                 aes(x=species, y= dispersion)) + 
  geom_point() +
  geom_line(group=F) + 
  #geom_bar(stat = "identity", aes(fill=species))  +
  guides(color = "none") +
  theme_bw() +  
  
  theme(axis.text.x =element_text(angle = 30))

write_csv(nfd.dispersion.df_sp,
          file.path("Results",
                    "nfd.dispersion.df_sp.csv"), 
          na = "NA", append = F,
          col_names =TRUE)

ggsave(plot = nfd.dispersion.plot_sp ,
       filename = paste0("nfd.dispersion.plot_sp",".pdf"),
       path = 'Figures',
       #scale=0.6,
       #width = 20, height = 10, dpi = 320, units = "in",
       device = "pdf")


library(RColorBrewer)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))


nfd.dispersion.plot <- ggplot(nfd.dispersion.df , 
                              aes(x=year, y= dispersion, 
                                  group= species, color= species)) + 
  geom_line(size=0.8, alpha=0.9) + 
  #geom_bar(stat = "identity", aes(fill=year))  +
  #guides(color = "none") +
  scale_color_manual(values=getPalette(length(unique(nfd.dispersion.df$species)))) +
  theme_bw() + 
  theme(axis.text.x =element_text(angle = 30))
ggplotly(nfd.dispersion.plot)

ggsave(plot = nfd.dispersion.plot ,
       filename = paste0("nfd.dispersion.plot",".pdf"),
       path = '/Users/lisabuche/Documents/Stage/Project/N-F-plane/Figures',
       #scale=0.6,
       #width = 20, height = 10, dpi = 320, units = "in",
       device = "pdf" )

write_csv(nfd.dispersion.df,
          file.path("Results",
                    "nfd.dispersion.df.csv"), 
          na = "NA", append = F,
          col_names =TRUE)


#######################################################################################################----
# to see the % of each Interaction 
#######################################################################################################----

Faciliation <- nrow(subset(NFD, NFD$ND >= 1)) / nrow(subset(NFD, !is.na(NFD$ND)))

Pos.freq.dep <- nrow(subset(NFD, NFD$ND < 0 )) / nrow(subset(NFD, !is.na(NFD$ND)))

Neg.freq.dep <- nrow(subset(NFD, NFD$ND <1 & NFD$ND >0 )) / nrow(subset(NFD, !is.na(NFD$ND)))

Neutrality <- nrow(subset(NFD, NFD$ND == 0 )) / nrow(subset(NFD, !is.na(NFD$ND)))

Faciliation + Pos.freq.dep + Neg.freq.dep + Neutrality

NFD.interaction <- data.frame( interaction = c("Faciliation", "Pos.freq.dep",
                                               "Neg.freq.dep","Neutrality"),
                               value = c(Faciliation, Pos.freq.dep,
                                         Neg.freq.dep,Neutrality))

NFD.interaction.plot <- ggplot(NFD.interaction, aes(x=interaction, y = value, fill= interaction)) +
  geom_bar(stat= "identity") + 
  theme_bw()
ggsave( plot= NFD.interaction.plot,
        filename = "NFD.interaction.pdf",
        path = "Figures/")

#######################################################################################################----
# Interactions outcome
#######################################################################################################----
NFD$outcome[i] <- NA
for ( i in c(1:nrow(NFD))){
  if (is.na(NFD$ND[i]) | is.na(NFD$FD[i])) next
  NFD$outcome[i] <- "Exclusion"
  if (NFD$ND[i] > NFD$FD[i]){
    NFD$outcome[i] <- "Coexistence"
  }
  
}

#========================================================================================================
#==========================
# Analyse of the dist in parallele to abundance and traits
#==========================
#========================================================================================================

plant_species_traits <- read.csv("data/plant_species_traits.csv")

precipitation_per_year <- read.csv("data/precipitation_per_year.csv")

nfd.dispersion.df_year <- spread(subset(nfd.dispersion.df, 
                                        select = c("species","year","dispersion")),
                                 key="year", "dispersion")
nfd.dispersion.df_nd_year <- spread(subset(nfd.dispersion.df, 
                                           select = c("species","year","nd")),
                                    key="year", nd)
nfd.dispersion.df_fd_year <- spread(subset(nfd.dispersion.df, 
                                           select = c("species","year","fd")),
                                    key="year", fd)


str(nfd.dispersion.df)
str(plant_species_traits)
names(plant_species_traits) <- c("species", names(plant_species_traits)[2:16])
str(precipitation_per_year)
names(precipitation_per_year) <- c("year", names(precipitation_per_year)[2:11])
precipitation_per_year$year <- as.character(precipitation_per_year$year)
nfd.dispersion.df_year_traits <- left_join(nfd.dispersion.df,plant_species_traits, by= c("species"))
nfd.dispersion.df_year_precipitation <- left_join(nfd.dispersion.df,precipitation_per_year, by= c("year"))
str(nfd.dispersion.df_year_traits)
nfd.dispersion.df_year_traits <- subset(nfd.dispersion.df_year_traits,
                                        select=c("dispersion",names(plant_species_traits)[5:16]))

nfd.dispersion.df_year_traits  <- scale(nfd.dispersion.df_year_traits , 
                                        center = TRUE, scale =TRUE) 
nfd.dispersion.df_year_traits  <- data.frame(nfd.dispersion.df_year_traits)
# Centrage (center = TRUE) et réduction (scale=TRUE) 
# des données du data.frame "mesures".
