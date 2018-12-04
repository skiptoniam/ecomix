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
fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli',
                   # weights = as.matrix(site_spp_weights),
                   n_mixtures = 4,
                   control = species_mix.control(minimum_sites_prevelance = 50,init_method = 'kmed'))

simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="poisson")
model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                covariate_data = simulated_data$covariate_data[,-1])
fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'poisson',
                   # weights = as.matrix(site_spp_weights),
                   n_mixtures = 4,
                   control = species_mix.control(minimum_sites_prevelance = 30,init_method = 'kmeans'))


simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="negative_binomial")
model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                covariate_data = simulated_data$covariate_data[,-1])
fm3 <- species_mix(sam_form, sp_form, model_data, distribution = 'negative_binomial',
                   # weights = as.matrix(site_spp_weights),
                   n_mixtures = 4,
                   control = species_mix.control(minimum_sites_prevelance = 50,init_method = 'kmed'))

simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="gaussian")
model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                covariate_data = simulated_data$covariate_data[,-1])
fm4 <- species_mix(sam_form, sp_form, model_data, distribution = 'gaussian',
                   # weights = as.matrix(site_spp_weights),
                   n_mixtures = 4,
                   control = species_mix.control(minimum_sites_prevelance = 50,init_method = 'kmed'))



simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="ippm")
model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                covariate_data = simulated_data$covariate_data[,-1])
simulated_data$background_weights -> wts
head(wts)
wts[!is.na(wts)&wts>0] <- (wts[!is.na(wts)&wts>0])*100

fm3 <- species_mix(sam_form, sp_form, model_data, distribution = 'ippm',
                   weights = wts,#simulated_data$background_weights,
                   n_mixtures = 4, standardise = FALSE,
                   control = species_mix.control(minimum_sites_prevelance = 50,init_method = 'kmeans'))
fm3
