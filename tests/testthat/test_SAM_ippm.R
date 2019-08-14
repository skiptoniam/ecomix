context('species_mix ippm')

library(ecomix)
library(raster)
library(scales)

testthat::test_that('species mix ippm', {

  # library(ecomix)
  # library(raster)
  rm(list=ls())
  set.seed(42)
  n_g <- 4
  n_sp <- 25
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:25),collapse = ','),")~1+x1+x2"))
  alphas <- runif(50,-8,-6)
  betas <- matrix(c(0.5,-1.3,
                     1.8,   1,
                     -1.2, 0.1),3,2,byrow=TRUE)

  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form, ~1, dat,n_mixtures = 3,
                                         alpha=alphas,beta = betas,
                                         distribution = "ippm")
  wts <- attr(simulated_data,"ippm_weights")
  y <- simulated_data[,1:n_sp]
  X <- simulated_data[,-1:-n_sp]
  W <- X[,1,drop=FALSE]
  X <- X[,-1]
  y_is_na <- is.na(y)
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- as.matrix(wts)
  G <- 3
  S <- n_sp
  ss <- 1
  disty <- 3
  inits <- NULL
  control <- species_mix.control(quiet = FALSE)
  offset <- rep(0,nrow(X))
  npx <- 2

  # test if one species ippm working - expect matrix of coefs back
  one_sp_ippm <- ecomix:::apply_glmnet_sam_inits(ss = ss, y = y, X = X, W = W,
                                        site_spp_weights = site_spp_weights,
                                        offset = offset, y_is_na = y_is_na,disty = disty)
  testthat::expect_is(one_sp_ippm,'list')

  # check that many species ippms work - expect back a list.
  all_sp_ippm <-surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits,
                                      y, X, W, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(all_sp_ippm,'list')

  alpha <- lapply(all_sp_ippm, `[[`, 1)
  testthat::expect_length(unlist(alpha),S)

  beta <- lapply(all_sp_ippm, `[[`, 2)
  testthat::expect_length(do.call(rbind, beta),S*npx)

  i <- 1
  testthat::expect_length(ecomix:::apply_glmnet_sam_inits(i, y, X, W, site_spp_weights, offset, y_is_na, disty),4)
  fm_ippmint <- surveillance::plapply(1:S, ecomix:::apply_glmnet_sam_inits,
                                                   y, X, W, site_spp_weights, offset,
                                                   y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_ippmint)[1,],S)

  #get the taus
  starting_values <- ecomix:::initiate_fit_sam(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
  fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,gamma=starting_values$gamma,theta=starting_values$theta)
  first_fit <- list(y = y, x = X, W = W, site_spp_weights=site_spp_weights,
                    offset=offset, y_is_na=y_is_na)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  taus <- ecomix:::shrink_taus(taus, max_tau=1/G + 0.1, G)

  # ## now let's try and fit the optimisation
  sp_form <- ~1
  fm1 <- species_mix(sam_form, sp_form, simulated_data,
                     distribution = 'ippm',
                     n_mixtures=3)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, simulated_data, distribution = 'ippm',
                     n_mixtures=3,control=species_mix.control(em_prefit = FALSE),
                     standardise = FALSE)
  testthat::expect_s3_class(fm2,'species_mix')


})
