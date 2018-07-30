context('species_mix poisson')

library(ecomix)


testthat::test_that('species mix poisson', {

set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="poisson")
model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                covariate_data = simulated_data$covariate_data[,-1])

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
fm_poissonint <- surveillance::plapply(1:S, ecomix:::apply_glmnet_poisson, y, X, weights, offset, .parallel = control$cores, .verbose = !control$quiet)
testthat::expect_length(do.call(cbind,fm_poissonint)[1,],S)

# test that the starting values work.
testthat::expect_length(tmp <- ecomix:::get_starting_values_poisson(y,X,offset,weights,G,S,control),10)

#get the taus
starting_values <- ecomix:::initiate_fit_poisson(y, X, weights, offset, G, S, control)
fits <- list(betas=starting_values$mix_coefs, alphas=starting_values$sp_intercepts)
first_fit <- list(x = X, y = y, weights=weights, offset=offset)

# get the loglikelihood based on these values
logls <- ecomix:::get_logls_poisson(first_fit, fits, G, S)
pis <- rep(1/G, G)
taus <- ecomix:::get_taus(pis, logls, G, S)
taus <- ecomix:::skrink_taus(taus, max_tau=1/G + 0.1, G)

## get to this in a bit
gg <- 1
testthat::expect_length(ecomix:::apply_glm_poisson_group_tau(gg, y, X, taus),2)

# ## now let's try and fit the optimisation
sv <- ecomix:::get_starting_values_poisson(y,X,offset,weights, G, S,control)
y_is_na <- is.na(y)
tmp <- ecomix:::sam_optimise(y,X,offset,sv$spp_wts,sv$site_spp_wts, y_is_na, sv$nS, sv$nG, sv$nObs, disty=2, start_vals = sv, control)
testthat::expect_length(tmp,15)

## most of the internal functions seem to be working.
## now let's test the species_mix function
sp_form <- ~1
fmp <- species_mix(sam_form, sp_form, model_data, distribution = 'poisson', n_mixtures=3, control = species_mix.control(quiet=TRUE,calculate_hessian_cpp = FALSE))
testthat::expect_s3_class(fmp, "species_mix")
testthat::expect_s3_class(fmp, "poisson")


})
