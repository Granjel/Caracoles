
# obtain interaction coefficients with a negative binomial model
# library(cxr)
library(tidyverse)
library(tidybayes)
library(brms)
library(rstan)
# source("./R/negbin_bayesian.R")

##################
# type of fit
# pool observations across space, time, both, or none

# fit.type <- "homogeneous" # pool together all plots and years
# fit.type <- "heterogeneous_space" # pool together different years, differentiate by plot
fit.type <- "heterogeneous_time" # pool together all plots, differentiate by year
# fit.type <- "heterogeneous_both" # independent fits for every plot and year
##################
# store results?
write.results <- T
##################

# read data
competition <- read.csv2(file = "data/competition.csv",header = T,stringsAsFactors = F)
sp.traits <- read.csv2(file = "data/plant_species_traits.csv",header = T,stringsAsFactors = F)

# subset species to those with independent germination estimates
valid.sp <- sp.traits$species.code[which(!is.na(sp.traits$germination.rate))]

# only positive number of seeds
competition <- subset(competition,seed > 0)

#focal species?
focal.sp <- sort(unique(competition$focal))

# keep only valid species that are also focal
# also, discard other neighbour species, so that we have square competition matrices
# the numbers of non-focal neighbours are negligible (see "neighbour_data_cheks.R")
all.sp <- focal.sp[focal.sp %in% valid.sp]

competition <- subset(competition, neighbour %in% all.sp)
competition <- subset(competition, focal %in% all.sp)

# spread the data from long to wide format
competition.data <- tidyr::spread(competition,key = neighbour,value = number,fill = 0)

# covariate: salinity
# salinity <- read.csv2(file = "../Caracoles/data/salinity.csv",header=T,stringsAsFactors = F)
# salinity$sum_salinity <- as.numeric(salinity$sum_salinity)
# salinity <- salinity[,c("plot","subplot","year","sum_salinity")]
# 
# # join the data
# full.data <- left_join(competition.data,salinity)
# 
# # which columns hold abundances of neighbours?
# neighbour.columns <- which(names(full.data) %in% sort(unique(competition$neighbour)))
# # how many neighbour species?
# num.neighbours <- length(neighbour.columns)
# # same for salinity
# covariate.columns <- which(names(full.data) == "sum_salinity")
# num.covariates <- 1

##############################
# initialize data structures
# in this case, a big dataframe in long format
# is the best option for analysis later on

all.ids <- sort(unique(competition$plot))
all.years <- sort(unique(competition$year))

interaction.coefs <- NULL
# interaction.coefs <- expand.grid(all.years,all.ids,all.sp,all.sp)
# names(interaction.coefs) <- c("year","plot","focal","neighbour")
# interaction.coefs$magnitude <- NA_real_

lambda.coefs <- NULL
# lambda.coefs <- expand.grid(all.years,all.ids,all.sp)
# names(lambda.coefs) <- c("year","plot","sp")
# lambda.coefs$lambda <- NA_real_

###############################
# main loop
# rescount <- 1
for(i.sp in 1:length(all.sp)){# each focal sp
  
  print("*********************************")
  print(paste(date()," - starting focal sp ",all.sp[i.sp],", model: ",fit.type,sep=""))
  
  sp.data <- subset(competition, focal == all.sp[i.sp])
  sp.data.long <- spread(sp.data,neighbour,number)

  my.data <- sp.data.long[,c("seed","year","plot",all.sp)]
  
  # which neighbours actually occur with this focal sp
  neighbours <- colnames(my.data)[!colnames(my.data) %in% c("seed","year","plot")]
  
  # remove species that are never observed to co-occur
  for(sp in names(which(colSums(my.data[,neighbours,drop=FALSE])==0))){
    my.data[,sp] <- NULL
  }
  # update list of neighbours
  neighbours <- colnames(my.data)[!colnames(my.data) %in% c("seed","year","plot")]
  
  # specify formula
  if(fit.type == "homogeneous"){
    model.formula <- paste0("seed ~ 1 + ",paste0(neighbours,collapse = " + "))
  }else if(fit.type == "heterogeneous_space"){
    model.formula <- paste0("seed ~ (1 + ",paste0(neighbours,collapse = " + "),"|plot)")
  }else if(fit.type == "heterogeneous_time"){
    model.formula <- paste0("seed ~ (1 + ",paste0(neighbours,collapse = " + "),"|year)")
  }else if(fit.type == "heterogeneous_both"){
    model.formula <- paste0("seed ~ (1 + ",paste0(neighbours,collapse = " + "),"|plot:year)")
  }

  # run model
  m <- try(brms::brm(as.formula(model.formula),
                     data = my.data, 
                     family = brms::negbinomial(),
                     iter = 5000,
                     chains = 6,
                     cores = 6,
                     control = list(adapt_delta = 0.99)
                     ))#,prior = my.prior,refresh = 0))
  
  # retrieve coefficients
  if(fit.type == "homogeneous"){
    lambda.draws <- m %>% spread_draws(b_Intercept)
    alpha.draws <- NULL
    
    lambda.draws <- lambda.draws[,3:4]
    names(lambda.draws) <- c("draw","log.lambda")
    
    # hack from https://github.com/mjskay/tidybayes/issues/38
    my.comp <- paste("b_",neighbours,sep="")
    alpha.draws <- m %>% spread_draws(!!!syms(my.comp))
    
    names(alpha.draws)[which(names(alpha.draws) %in% my.comp)] <- neighbours
    alpha.draws <- tidyr::gather(alpha.draws,neighbour,magnitude,-.chain,-.iteration,-.draw)
    alpha.draws <- alpha.draws[,c(".draw","neighbour","magnitude")]
    names(alpha.draws)[1] <- "draw"
    
    ##### NOTE
    # the sign of alpha coefficients is opposite to that of the statistical model
    # (see code from mayfield & stouffer 2017)
    # STORE ALL RAW DRAWS WITH SIGNS ALREADY SWITCHED
    alpha.draws$magnitude <- -alpha.draws$magnitude
    ######
    
    # this is the same as averaging over draws
    my.coefs <- brms::fixef(m)
    lambda <- my.coefs["Intercept","Estimate"]
    
    ###### SEE NOTE ABOVE
    alphas <- - my.coefs[neighbours,"Estimate"]
    ######
    
    mean.lambdas <- expand.grid(year = all.years,plot = all.ids)
    mean.lambdas$log.lambda <- lambda
    
    mean.alphas <- expand.grid(year = all.years,plot = all.ids,neighbour = neighbours)
    mean.alphas$mean.magnitude <- alphas[match(mean.alphas$neighbour,names(alphas))]
    
    # add current focal sp
    mean.lambdas$focal <- all.sp[i.sp]
    mean.alphas$focal <- all.sp[i.sp]
    
  }else if(fit.type == "heterogeneous_space"){
    
    all.draws <- m %>% spread_draws(b_Intercept,r_plot[plot,term])
    lambda.draws <- all.draws[all.draws$term == "Intercept",c(3,5,4,7)]
    names(lambda.draws) <- c("draw","plot","log.lambda","log.lambda_offset")
    mean.lambdas <- lambda.draws %>% group_by(plot) %>% 
      summarise(mean.log.lambda = mean(log.lambda), 
                mean.log.offset = mean(log.lambda_offset))
    
    alpha.draws <- all.draws[all.draws$term != "Intercept",c(3,5,6,7)]
    names(alpha.draws) <- c("draw","plot","neighbour","magnitude")
    
    ##### NOTE
    # the sign of alpha coefficients is opposite to that of the statistical model
    # (see code from mayfield & stouffer 2017)
    # STORE ALL RAW DRAWS WITH SIGNS ALREADY SWITCHED
    alpha.draws$magnitude <- -alpha.draws$magnitude
    ######
    
    mean.alphas <- alpha.draws %>% group_by(plot,neighbour) %>% 
      summarise(mean.magnitude = mean(magnitude))
    
    # add year field
    mean.lambdas <- mean.lambdas[,c("plot","mean.log.lambda","mean.log.offset")]
    all.plots.lambdas <- left_join(unique(my.data[,c("year","plot")]),mean.lambdas)
    
    # name consistency
    mean.lambdas <- dplyr::arrange(all.plots.lambdas,year,plot)
    
    # same as for lambda, add year field
    all.plots.alphas <- left_join(unique(my.data[,c("year","plot")]),mean.alphas)
    my.data.long <- gather(my.data,key = "neighbour",value = "number",-seed,-year,-plot)
    all.plots.comp <- left_join(my.data.long,mean.alphas)
    all.plots.comp <- unique(all.plots.comp[,c("year","plot","neighbour","number","mean.magnitude")])
    all.plots.comp <- subset(all.plots.comp,number > 0)
    all.plots.alphas <- unique(all.plots.comp[,c("year","plot","neighbour","mean.magnitude")])
    
    # name consistency
    mean.alphas <- dplyr::arrange(all.plots.alphas,year,plot,neighbour)
    
    # add current focal sp
    mean.lambdas$focal <- all.sp[i.sp]
    mean.alphas$focal <- all.sp[i.sp]
    
  }else if(fit.type == "heterogeneous_time"){
    
    all.draws <- m %>% spread_draws(b_Intercept,r_year[year,term])
    lambda.draws <- all.draws[all.draws$term == "Intercept",c(3,5,4,7)]
    names(lambda.draws) <- c("draw","year","log.lambda","log.lambda_offset")
    mean.lambdas <- lambda.draws %>% group_by(year) %>% 
      summarise(mean.log.lambda = mean(log.lambda), 
                mean.log.offset = mean(log.lambda_offset))
    
    alpha.draws <- all.draws[all.draws$term != "Intercept",c(3,5,6,7)]
    names(alpha.draws) <- c("draw","year","neighbour","magnitude")
    
    ##### NOTE
    # the sign of alpha coefficients is opposite to that of the statistical model
    # (see code from mayfield & stouffer 2017)
    # STORE ALL RAW DRAWS WITH SIGNS ALREADY SWITCHED
    alpha.draws$magnitude <- -alpha.draws$magnitude
    ######
    
    mean.alphas <- alpha.draws %>% group_by(year,neighbour) %>% 
      summarise(mean.magnitude = mean(magnitude))

    # add plot field
    mean.lambdas <- mean.lambdas[,c("year","mean.log.lambda","mean.log.offset")]
    all.plots.lambdas <- left_join(unique(my.data[,c("year","plot")]),mean.lambdas)
    
    # name consistency
    mean.lambdas <- dplyr::arrange(all.plots.lambdas,year,plot)
    
    # same as for lambda, add plot field
    
    all.plots.alphas <- left_join(unique(my.data[,c("year","plot")]),mean.alphas)
    my.data.long <- gather(my.data,key = "neighbour",value = "number",-seed,-year,-plot)
    all.plots.comp <- left_join(my.data.long,mean.alphas)
    all.plots.comp <- unique(all.plots.comp[,c("year","plot","neighbour","number","mean.magnitude")])
    all.plots.comp <- subset(all.plots.comp,number > 0)
    all.plots.alphas <- unique(all.plots.comp[,c("year","plot","neighbour","mean.magnitude")])
    
    # name consistency
    mean.alphas <- dplyr::arrange(all.plots.alphas,year,plot,neighbour)
    
    # add current focal sp
    mean.lambdas$focal <- all.sp[i.sp]
    mean.alphas$focal <- all.sp[i.sp]
    
  }else if(fit.type == "heterogeneous_both"){
    all.draws <- m %>% spread_draws(b_Intercept,`r_plot:year`[`plot:year`,term])
    lambda.draws <- all.draws[all.draws$term == "Intercept",c(3,5,4,7)]
    names(lambda.draws) <- c("draw","site","log.lambda","log.lambda_offset")
    mean.lambdas <- lambda.draws %>% group_by(site) %>% 
      summarise(mean.log.lambda = mean(log.lambda), 
                mean.log.offset = mean(log.lambda_offset))
    
    alpha.draws <- all.draws[all.draws$term != "Intercept",c(3,5,6,7)]
    names(alpha.draws) <- c("draw","site","neighbour","magnitude")
    
    ##### NOTE
    # the sign of alpha coefficients is opposite to that of the statistical model
    # (see code from mayfield & stouffer 2017)
    # STORE ALL RAW DRAWS WITH SIGNS ALREADY SWITCHED
    alpha.draws$magnitude <- -alpha.draws$magnitude
    ######
    
    mean.alphas <- alpha.draws %>% group_by(site,neighbour) %>% summarise(mean.magnitude = mean(magnitude))
    
    mean.lambdas$plot <- substr(mean.lambdas$site,1,1)
    mean.lambdas$year <- substr(mean.lambdas$site,3,6)
    mean.lambdas <- mean.lambdas[,c("year","plot","mean.log.lambda","mean.log.offset")]
    
    mean.alphas$plot <- substr(mean.alphas$site,1,1)
    mean.alphas$year <- substr(mean.alphas$site,3,6)
    mean.alphas <- mean.alphas[,c("year","plot","neighbour","mean.magnitude")]
    
    # add current focal sp
    mean.lambdas$focal <- all.sp[i.sp]
    mean.alphas$focal <- all.sp[i.sp]
    
  }
  
  # append to results
  lambda.coefs <- rbind(lambda.coefs,mean.lambdas)
  interaction.coefs <- rbind(interaction.coefs,mean.alphas)
  
  # store species posterior samples
  write.csv2(lambda.draws,paste("./results/temp/lambda_samples_",fit.type,"_",all.sp[i.sp],".csv",sep=""),row.names = F)
  write.csv2(alpha.draws,paste("./results/temp/alpha_samples_",fit.type,"_",all.sp[i.sp],".csv",sep=""),row.names = F)

  # store species means
  write.csv2(mean.lambdas,paste("./results/temp/lambda_coefs_",fit.type,"_",all.sp[i.sp],".csv",sep=""),row.names = F)
  write.csv2(mean.alphas,paste("./results/temp/alpha_coefs_",fit.type,"_",all.sp[i.sp],".csv",sep=""),row.names = F)
  
}# for each focal sp

# get back model parameters from regression coefficients
lambda.coefs$year <- as.numeric(lambda.coefs$year)
lambda.coefs$plot <- as.numeric(lambda.coefs$plot)
if(fit.type != "homogeneous"){
  lambda.coefs$lambda <- exp(lambda.coefs$mean.log.lambda + lambda.coefs$mean.log.offset)
}else{
  lambda.coefs$lambda <- exp(lambda.coefs$log.lambda)
}
lambda.coefs <- lambda.coefs[,c("year","plot","focal","lambda")]
names(lambda.coefs)[3] <- "sp"


# SIGN OF ALPHA COEFS -----------------------------------------------------
# alphas are ok (signs have been switched above)
# to summarise, as shown by Stouffer et al.:
# alphas in the negative binomial model are: competition > 0, facilitation < 0
# alphas returned by the regression model are: competition < 0, facilitation > 0
# so, switch them. Actually, then, the alphas stored are the inverse of the intuitive

# fill up remaining positions with NA
all.positions <- expand.grid(year = all.years,plot = all.ids,sp = all.sp)
all.positions$lambda.na <- NA_real_
all.positions <- left_join(all.positions,lambda.coefs)
all.lambda <- all.positions[,c("year","plot","sp","lambda")]

interaction.coefs$year <- as.numeric(interaction.coefs$year)
interaction.coefs$plot <- as.numeric(interaction.coefs$plot)
alpha.positions <- expand.grid(year = all.years,plot = all.ids,focal = all.sp,neighbour = all.sp)
alpha.positions$m <- NA_real_
all.alpha <- left_join(alpha.positions,interaction.coefs)

all.alpha <- all.alpha[,c("year","plot","focal","neighbour","mean.magnitude")]
names(all.alpha)[5] <- "magnitude"

# # 
if(write.results){
      write.csv2(all.lambda,file = paste("./results/lambda_negbin_multilevel_",fit.type,"_confirmatory.csv",sep=""),row.names = F)
      write.csv2(all.alpha,file = paste("./results/alpha_negbin_multilevel_",fit.type,"_confirmatory.csv",sep=""),row.names = F)
}

# small test to see if models return proper estimates
# -> They do :)

# lambda.coefs.orig <- read.csv2(file = "./results/lambda_negbin_heterogeneous_both.csv",header=T)
# comp.long <- spread(competition,neighbour,number)
# 
# for(i.sp in 1:length(all.sp)){# each focal sp
  # for(i.year in 1:length(all.years)){
  #   for(i.plot in 1:length(all.ids)){
  #     my.num <- sum(comp.long$focal == all.sp[i.sp] &
  #                     comp.long$year == all.years[i.year] &
  #                     comp.long$plot == all.ids[i.plot])
  #     my.lambda <- lambda.coefs.orig$lambda[lambda.coefs$year == all.years[i.year] &
  #                                        lambda.coefs$plot == all.ids[i.plot] &
  #                                        lambda.coefs$sp == all.sp[i.sp]]
  #     message(paste("sp ",all.sp[i.sp],", plot ",all.ids[i.plot],", year ",all.years[i.year],", with n=",my.num," has lambda=",round(my.lambda,2),sep = ""))
  #     # if(is.na(my.lambda) & my.num < 5){
  #     #   message(paste("sp ",all.sp[i.sp],", plot ",all.ids[i.plot],", year ",all.years[i.year],", with n=",my.num," not enough data",sep=""))
  #     # }
  #     # if(is.na(my.lambda) & my.num > 4){
  #     #   message(paste("sp ",all.sp[i.sp],", plot ",all.ids[i.plot],", year ",all.years[i.year],", with n=",my.num," not calculated",sep=""))
  #     # }
  #   }
  # }
# }

