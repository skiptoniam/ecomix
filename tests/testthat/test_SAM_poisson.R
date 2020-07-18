context('species_mix generic functions five: poisson functions')
library(ecomix)

testthat::test_that('species mix poisson', {

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~x1+x2"))
  alpha <- rnorm(20,1, 0.5)
  beta <- matrix(c(-2.6,0.5,
                   -0.9,1.0,
                   0.9,-1.9),
                 3,2,byrow=TRUE)
  dat <- data.frame(y=1, x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  simulated_data <- species_mix.simulate(archetype_formula = sam_form,
                                         species_formula = ~1,
                                         dat = dat,
                                         n_mixtures = 3,
                                         alpha=alpha,
                                         beta=beta,
                                         distribution = "poisson")
  y <- as.matrix(simulated_data[,grep("spp",colnames(simulated_data))])
  X <- simulated_data[,-grep("spp",colnames(simulated_data))]
  W <- as.matrix(X[,1,drop=FALSE])
  X <- as.matrix(X[,-1])
  offset <- rep(0,nrow(y))
  weights <- rep(1,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- 3
  S <- ncol(y)
  control <- species_mix.control()
  disty <- 2
  size <- rep(1,nrow(y))

  # test a single gaussian model
  i <- 1
  testthat::expect_length(ecomix:::apply_glmnet_sam_inits(i, y, X, W, site_spp_weights, offset, y_is_na, disty,size),4)
  fm_poissonint <- surveillance::plapply(1:S, ecomix:::apply_glmnet_sam_inits,
                                         y, X, W, site_spp_weights, offset,
                                         y_is_na, disty, size, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_poissonint)[1,],S)

  #get the taus
  starting_values <- ecomix:::initiate_fit_sam(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, size, control)
  fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,gamma=starting_values$gamma,theta=starting_values$theta)
  first_fit <- list(y = y, x = X, W = W, weights=site_spp_weights, offset=offset,size=size)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  taus <- ecomix:::shrink_taus(taus, max_tau=1/G + 0.1, G)

  ## get to this in a bit
  gg <- 1
  testthat::expect_length(ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W, site_spp_weights, offset, y_is_na, disty, taus, fits, logls$fitted, size),2)

  # ## now let's try and fit the optimisation
  sv <- ecomix:::get_starting_values_sam(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, size, control)
  tmp <- ecomix:::sam_optimise(y, X, W, offset, spp_weights, site_spp_weights, y_is_na, S, G, disty, size, start_vals = sv, control)
  testthat::expect_length(tmp,18)

  sp_form <- ~1
  fm1 <- species_mix(sam_form, sp_form, simulated_data, distribution = 'poisson',
                     n_mixtures=3)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, simulated_data, distribution = 'poisson',
                     n_mixtures=3,control=species_mix.control(em_prefit = FALSE),
                     standardise = FALSE)
  testthat::expect_s3_class(fm2,'species_mix')
})

