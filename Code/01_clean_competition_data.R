
# read and clean data

comp <- read.csv("/home/david/Work/Projects/EBD/Caracoles/data/competition.csv")

# remove some species with few observations as focal (also as neighbours)
# with Melilotus, do we join them all together or keep them separate?

to.remove <- c("ANAR","ACHI","ARTE","NEWGRASS","COSQ","CRCR","DA","LYTR","MEPO",
               "RAPE")

# -------------------------------------------------------------------------

comp.clean <- subset(comp, !focal %in% to.remove)
comp.clean[,to.remove] <- NULL

focals <- unique(comp.clean$focal)
neighs <- names(comp.clean)[8:26]

# FRPU is present as neighbour but not focal, 
# and has only 27 observations as neighbour

comp.clean$FRPU <- NULL

# now, the list of focals and neighbours should be consistent

focals <- unique(comp.clean$focal)
neighs <- names(comp.clean)[8:25]

comp.clean$seed <- round(comp.clean$seed)

write.csv(comp.clean,"/home/david/Work/Projects/EBD/population_niche/data/competition_clean.csv")

