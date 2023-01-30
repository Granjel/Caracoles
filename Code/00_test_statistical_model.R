library(tidyverse)
library(rbms)
library(tidybayes)

# simulate species effects on individual fecundity

lambda.df <- expand.grid(sp = c("sp1","sp2","sp3"), 
                         site = c("a","b"), lambda = 0)

lambda.df$lambda[lambda.df$sp == "sp1"] <- 50
lambda.df$lambda[lambda.df$sp == "sp2" & lambda.df$site == "a"] <- 100
lambda.df$lambda[lambda.df$sp == "sp2" & lambda.df$site == "b"] <- 30
lambda.df$lambda[lambda.df$sp == "sp3" & lambda.df$site == "a"] <- 10
lambda.df$lambda[lambda.df$sp == "sp3" & lambda.df$site == "b"] <- 40

a1 <- matrix(data = c(.01,.07,.02,.01,.08,.01,.03,0,.04),nrow = 3,byrow = TRUE)
a2 <- matrix(data = c(.01,.03,.03,.02,.04,.04,.03,.01,.06),nrow = 3,byrow = TRUE)

site1 <- data.frame(focal.sp = rep(c("sp1","sp2","sp3"),100), 
                     seed = 0, 
                     sp1 = sample(1:20,300,replace = TRUE),
                     sp2 = sample(1:100,300,replace = TRUE),
                     sp3 = sample(1:30,300,replace = TRUE),
                    site = "a")
site2 <- data.frame(focal.sp = rep(c("sp1","sp2","sp3"),100), 
                    seed = 0, 
                    sp1 = sample(1:10,300,replace = TRUE),
                    sp2 = sample(1:10,300,replace = TRUE),
                    sp3 = sample(1:50,300,replace = TRUE),
                    site = "b")

# site1$seed <- sapply(site1$focal.sp,FUN = function(x) lambda.df$lambda[lambda.df$sp == x] * exp(-()))

for(i.obs in 1:nrow(site1)){
  my.sp <- site1$focal.sp[i.obs]
  my.sp.num <- as.numeric(substr(my.sp,3,3))
  my.coefs <- sum(site1[i.obs,3:5] * a1[my.sp.num,])
  
  site1$seed[i.obs] <- (lambda.df$lambda[lambda.df$sp == my.sp &
                                          lambda.df$site == "a"] * exp(-my.coefs)) + rnorm(1,0,.2)
}

for(i.obs in 1:nrow(site2)){
  my.sp <- site2$focal.sp[i.obs]
  my.sp.num <- as.numeric(substr(my.sp,3,3))
  my.coefs <- sum(site2[i.obs,3:5] * a2[my.sp.num,])
  
  site2$seed[i.obs] <- (lambda.df$lambda[lambda.df$sp == my.sp &
                                           lambda.df$site == "b"] * exp(-my.coefs)) + rnorm(1,0,2)
}

site1$seed[site1$seed < 0] <- 0
site1$seed <- round(site1$seed)
site2$seed <- round(site2$seed)

hist(site1$seed)
hist(site2$seed)

sim.data <- rbind(site1,site2)

# -------------------------------------------------------------------------

focal <- "sp1"
my.data <- subset(sim.data,focal.sp == focal)
neighbours <- c("sp1","sp2","sp3")

m1.1 <-  paste0("seed ~ (1 + ",paste0(neighbours,collapse = " + "),"|site)")
# m1.2 <- "seed ~ sp1 + sp2 + sp3 + (1|site)"

m1 <- try(brms::brm(as.formula(m1.2),
                    data = my.data, 
                    family = brms::negbinomial(),
                    iter = 5000,
                    chains = 6,
                    cores = 6,
                    control = list(adapt_delta = 0.99)
))#,prior = my.prior,refresh = 0))

all.draws <- m1 %>% spread_draws(b_Intercept,r_site[site,term])
lambda.draws <- all.draws[all.draws$term == "Intercept",c(3,5,4,7)]
names(lambda.draws) <- c("draw","site","log.lambda","log.lambda_offset")
mean.lambdas <- lambda.draws %>% group_by(site) %>% 
  summarise(mean.log.lambda = mean(log.lambda), 
            mean.log.offset = mean(log.lambda_offset))

alpha.draws <- all.draws[all.draws$term != "Intercept",c(3,5,6,7)]
names(alpha.draws) <- c("draw","site","neighbour","magnitude")

alpha.draws$magnitude <- -alpha.draws$magnitude

mean.alphas <- alpha.draws %>% group_by(site,neighbour) %>% 
  summarise(mean.magnitude = mean(magnitude))

# -------------------------------------------------------------------------

# real lambda
lambda.df$lambda[lambda.df$sp == focal]
# expected lambda
exp(mean.lambdas$mean.log.lambda + mean.lambdas$mean.log.offset)

# real alphas, site 1
a1[1,]
# expected alphas, site 1
subset(mean.alphas, site == "a")

# real alphas, site 2
a2[1,]
# expected alphas, site 2
subset(mean.alphas, site == "b")
