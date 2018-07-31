context('species_mix bernoulli')
library(ecomix)


testthat::test_that('species mix bernoulli', {


  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2+z1+z2"))
  sp_form <- NULL
  theta <- matrix(c(1,-2.9,3.6,.5,-1,1,-0.9,1,-3,2,1,-.9,7.9,3,-2),3,5,byrow=TRUE)
  x1<-runif(100,0,2.5)
  z1<-rnorm(100,0,2.5)
  dat <- data.frame(y=rep(1,100),x1,x2=x1^2,z1,z2=z1^2)
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
                                               dat,theta,dist="bernoulli")
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

  # test a single bernoulli model
  ss <- 1
  testthat::expect_length(ecomix:::apply_glmnet_bernoulli_sp(ss, y, X, weights, offset),5)
  fm_bern_int <- surveillance::plapply(1:S, ecomix:::apply_glmnet_bernoulli_sp, y, X, weights, offset, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_bern_int)[1,],S)

  # test that the starting values work.
  testthat::expect_length(tmp <- ecomix:::get_starting_values_bernoulli_sp(y, X, weights, offset, G, S, control),11)

  #get the taus
  starting_values <- ecomix:::initiate_fit_bernoulli_sp(y, X, weights, offset, G, S, control)
  fits <- list(betas=starting_values$mix_coefs, alphas=starting_values$sp_intercepts)
  first_fit <- list(x = X, y = y, weights=weights, offset=offset)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_bernoulli_sp(first_fit, fits, G, S)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls, G, S)
  taus <- ecomix:::skrink_taus(taus, max_tau=1/G + 0.1, G)

  ## get to this in a bit
  gg <- 1
  testthat::expect_length(ecomix:::apply_glm_bernoulli_group_tau(gg, y, X, taus),4)

  # ## now let's try and fit the optimisation
  sv <- ecomix:::get_starting_values_bernoulli_sp(y, X, weights, offset, G, S, control)
  y_is_na <- is.na(y)
  tmp <- ecomix:::sam_optimise(y,X,offset,sv$spp_wts,sv$site_spp_wts, y_is_na, sv$nS, sv$nG, sv$nObs, disty=1, start_vals = sv, control)
  testthat::expect_length(tmp,15)

  ## most of the internal functions seem to be working.
  ## now let's test the species_mix function
  sp_form <- ~1
  fmb <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli', n_mixtures=3, control = species_mix.control(quiet=FALSE,calculate_hessian_cpp = FALSE))
  testthat::expect_s3_class(fmp, "species_mix")
  testthat::expect_s3_class(fmp, "bernoulli")

})
