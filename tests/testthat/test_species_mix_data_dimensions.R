context('species_mix data dimensions')

library(ecomix)
library(raster)
library(scales)

testthat::test_that('species mix one covariate and one group', {

set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
theta <- matrix(c(1,-1.6,-2.5,
                  1,-0.9,1,
                  1,2.9,2.9,
                  1,0.2,-0.4),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="bernoulli")
model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                covariate_data = simulated_data$covariate_data[,-1])

sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1"))
sp_form <- ~ 1
fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli',
                   # weights = as.matrix(site_spp_weights),
                   n_mixtures = 4,
                   control = species_mix.control(minimum_sites_prevelance = 50,init_method = 'kmeans'))


# test with one group
set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
theta <- matrix(c(1,-1.6,-2.5,
                  1,-0.9,1,
                  1,2.9,2.9,
                  1,0.2,-0.4),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="bernoulli")
model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                covariate_data = simulated_data$covariate_data[,-1])

sp_form <- ~ 1
fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli',
                   # weights = as.matrix(site_spp_weights),
                   n_mixtures = 1,
                   control = species_mix.control(minimum_sites_prevelance = 50,init_method = 'kmed'))

})
