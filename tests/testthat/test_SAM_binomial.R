context('species_mix generic functions two: binomial functions')
library(ecomix)

testthat::test_that('species mix bernoulii functions work', {

  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),
                                                      collapse = ','),
                                       ")~1+x1+x2"))
  sp_form <- ~ 1
  beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,7.9),3,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),
                    x2=stats::rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  model_data <- species_mix.simulate(archetype_formula=sam_form,
                                     species_formula=sp_form,
                                     dat,beta=beta,
                                     size = rep(100,nrow(dat)),
                                     distribution = "binomial")
  testthat::expect_message(fm1 <- species_mix(NULL, sp_form,
                                              model_data,
                                              distribution = 'binomial',
                                              n_mixtures=3))

  dup_spp_data <- cbind('spp1'=model_data[,1],model_data)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',c(1,1:20)),collapse = ','),")~1+x1+x2"))

  testthat::expect_message(fm1 <- species_mix(sam_form, sp_form, dup_spp_data,
                                              distribution = 'binomial',
                                              n_mixtures=3))
})


testthat::test_that('species mix binomial', {


  rm(list = ls())
  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),
                                                      collapse = ','),
                                       ")~1+x1+x2"))
  sp_form <- ~ 1
  beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,7.9),3,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),
                    x2=stats::rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(archetype_formula=sam_form,
                                     species_formula=sp_form,
                                     dat,beta=beta,
                                     size = rep(100,nrow(dat)),
                                     distribution = "binomial")
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
  disty <- 7
  size <- attr(simulated_data,"size")


  # test a single binomial model
  i <- 1
  testthat::expect_length(ecomix:::apply_glm_sam_inits(i, y, X, W, site_spp_weights, offset, y_is_na, disty, size),4)
  fm_binomialint <- surveillance::plapply(1:S, ecomix:::apply_glm_sam_inits,
                                           y, X, W, site_spp_weights, offset,
                                           y_is_na, disty, size, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_binomialint)[1,],S)

  #get the taus
  starting_values <- ecomix:::initiate_fit_sam(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, size, control)
  fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,gamma=starting_values$gamma,theta=starting_values$theta)
  first_fit <- list(y = y, x = X, W = W, weights=site_spp_weights, offset=offset, size=size)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  taus <- ecomix:::shrink_taus(taus, max_tau=1/G + 0.1, G)

  ## get to this in a bit
  gg <- 1
  testthat::expect_length(ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W, site_spp_weights, offset, y_is_na, disty, taus, fits, logls$fitted, size),2)

  # ## now let's try and fit the optimisation
  sv <- ecomix:::get_starting_values_sam(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, size, control=species_mix.control(em_steps = 100))
  tmp <- ecomix:::sam_optimise(y, X, W, offset, spp_weights, site_spp_weights, y_is_na, S, G, disty, size, start_vals = sv, control = species_mix.control())
  testthat::expect_length(tmp,18)

  ## test species mix fit
  inits <- NULL
  # tmp <- ecomix:::species_mix.fit(y=y, X=X, W=W, G=G, S=S,
  #                        spp_weights=spp_weights,
  #                        site_spp_weights=site_spp_weights,
  #                        offset=offset, disty=disty, y_is_na=y_is_na,
  #                        control=control, inits=inits)
  set.seed(123)
  tmp1 <- ecomix:::get_starting_values_sam(y, X, W, spp_weights, site_spp_weights,
                                           offset, y_is_na, G, S, disty,size,
                                           control=species_mix.control(em_refit = 5, em_steps = 10))
  set.seed(123)
  tmp2 <- ecomix:::fitmix_ECM_sam(y=y, X=X, W=W, G=G, S=S,
                         spp_weights=spp_weights,
                         site_spp_weights=site_spp_weights,
                         offset=offset, disty=disty, y_is_na=y_is_na,size = size,
                         control=species_mix.control(em_refit = 5, em_steps = 5))
  set.seed(123)
  tmp3 <- ecomix:::get_starting_values_sam(y, X, W, spp_weights, site_spp_weights,
                                           offset, y_is_na, G, S, disty,size,
                                           control = species_mix.control(em_prefit = FALSE))
  set.seed(123)
  tmp <- ecomix:::species_mix.fit(y=y, X=X, W=W, G=G, S=S,
                         spp_weights=spp_weights,
                         site_spp_weights=site_spp_weights,
                         offset=offset, disty=disty,size = size,y_is_na=y_is_na,
                         control=species_mix.control(print_cpp_start_vals = TRUE), inits=inits)

  sp_form <- ~1
  fm1 <- species_mix(sam_form, sp_form, simulated_data, distribution = 'binomial',
                       n_mixtures=3, size = size)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, simulated_data, distribution = 'binomial',
                     n_mixtures=3, size = size, control=species_mix.control(em_prefit = FALSE),
                     standardise = FALSE)
  testthat::expect_s3_class(fm2,'species_mix')

})

