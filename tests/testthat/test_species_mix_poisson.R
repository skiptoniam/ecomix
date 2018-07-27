context('species_mix poisson')

library(ecomix)


testthat::test_that('species mix poisson', {

set.seed(42)
form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- simulate_species_mix_data(form,~1,dat,theta,dist="poisson")

y <- simulated_data$species_data
X <- simulated_data$covariate_data
offset <- rep(0,nrow(y))
weights <- rep(1,nrow(y))
spp_wts <- rep(1,ncol(y))
site_spp_wts <- matrix(1,nrow(y),ncol(y))
y_is_na <- matrix(FALSE,nrow(y),ncol(y))
G <- length(simulated_data$pi)
S <- length(simulated_data$sp.int)
control <- species_mix.control()

# test a single poisson model
i <- 1
testthat::expect_length(ecomix:::apply_glmnet_poisson(i, y, X, weights, offset),3)

# test that the starting values work.
testthat::expect_length(tmp <- get_starting_values_poisson(y,X,offset,weights,S,G,control),10)

#get the taus
starting_values <- ecomix:::initiate_fit_poisson(y, X, weights, offset, G, S, control)
fits <- list(betas=starting_values$mix_coefs, alphas=starting_values$sp_intercepts)
first_fit <- list(x = X, y = y, weights=weights, offset=offset)

# get the loglikelihood based on these values
logls <- ecomix:::get_logls_poisson(first_fit, fits, G, S)
pis <- rep(1/G, G)
taus <- get_taus(pis, logls, G, S)
taus <- skrink_taus(taus, max_tau=1/G + 0.1, G)

## get to this in a bit
gg <- 1
testthat::expect_length(ecomix:::apply_glm_poisson_group_tau(gg, y, X, taus),2)

# ## now let's try and fit the optimisation
tmp <- get_starting_values_poisson(y,X,offset,weights,S,G,control)
y_is_na <- is.na(y)
res <- ecomix:::sam_optimise(y,X,offset,tmp$spp_wts,tmp$site_spp_wts, y_is_na, tmp$nS, tmp$nG, tmp$nObs, disty=2, start_vals = tmp, control)
testthat::expect_length(res,14)
})
