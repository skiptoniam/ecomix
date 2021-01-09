context('species_mix generic functions two: bernoulli functions')
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
                                     species_formula=sp_form,data = dat,
                                     beta=beta,family="bernoulli")
  testthat::expect_message(fm1 <- species_mix(NULL, sp_form,data = model_data,
                                              family = 'bernoulli',
                                              nArchetypes = 3))

  dup_spp_data <- cbind('spp1'=model_data[,1],model_data)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',c(1,1:20)),collapse = ','),")~1+x1+x2"))

  testthat::expect_message(fm1 <- species_mix(sam_form, sp_form, data=dup_spp_data,
                                              family = 'bernoulli',
                                              nArchetypes = 3))
})


testthat::test_that('species mix bernoulli', {

 library(ecomix)
  rm(list = ls())
  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:50),
                                                      collapse = ','),
                                       ")~1+x1+x2"))
  sp_form <- ~ 1
  beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,7.9),3,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,200),x1=stats::runif(200,0,2.5),
                    x2=stats::rnorm(200,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(archetype_formula=sam_form,
                                     species_formula=sp_form,data = dat,
                                     beta=beta,family="bernoulli")
  y <- as.matrix(simulated_data[,grep("spp",colnames(simulated_data))])
  X <- simulated_data[,-grep("spp",colnames(simulated_data))]
  W <- as.matrix(X[,1,drop=FALSE])
  X <- as.matrix(X[,-1])
  U <- NULL
  offset <- rep(0,nrow(y))
  weights <- rep(1,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- 3
  S <- ncol(y)
  control <- ecomix:::set_control_sam(list())#species_mix.control()
  disty <- 1
  size <- rep(1,nrow(y))
  powers <- rep(1.5,S)#attr(simulated_data,"powers") # yeah baby


  # test a single bernoulli model
  i <- 1
  testthat::expect_length(ecomix:::apply_species_fits(i, y, X, W, U, site_spp_weights, offset, y_is_na, disty, size, powers),5)
  fm_bernoulliint <- ecomix:::plapply(1:S, ecomix:::apply_species_fits,
                                           y, X, W, U, site_spp_weights, offset, y_is_na, disty, size, powers, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_bernoulliint)[1,],S)

  #get the taus
  sv <- ecomix:::get_initial_values_sam(y, as.data.frame(X), as.data.frame(W), U, site_spp_weights, offset, y_is_na, G, S, disty, size, powers, control)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(y, X, W, U, G, S, spp_weights,
                                  site_spp_weights, offset, y_is_na, disty,
                                  size, powers, control, sv, get_fitted = FALSE)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  taus <- ecomix:::shrink_taus(taus, G)

  ## get to this in a bit
  # gg <- 1
  # testthat::expect_length(ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W, site_spp_weights, offset, y_is_na, disty, taus, fits, logls$fitted, size),2)

  # ## now let's try and fit the optimisation
  start_vals <- ecomix:::starting_values_wrapper(y, as.data.frame(X), as.data.frame(W), U, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, size, powers, control)
  tmp <- ecomix:::sam_optimise(y, X, W, U, offset, spp_weights, site_spp_weights, y_is_na, S, G, disty, size, powers, start_vals = start_vals, control)
  testthat::expect_length(tmp,21)

  set.seed(123)
  tmp <- ecomix:::species_mix.fit(y=y, X=as.data.frame(X), W=as.data.frame(W), U=U, G=G, S=S,
                         spp_weights=spp_weights,
                         site_spp_weights=site_spp_weights,
                         offset=offset, disty=disty, y_is_na=y_is_na, size=size, powers=powers,
                         control=ecomix:::set_control_sam(list(print_cpp_start_vals = TRUE)))

  sp_form <- ~1
  fm1 <- species_mix(archetype_formula = sam_form, species_formula = sp_form,
                     data = simulated_data, family = 'bernoulli',
                     nArchetypes = 3)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, data = simulated_data, family = 'bernoulli',
                     nArchetypes = 3,control=list(em_prefit = FALSE),
                     standardise = FALSE)
  testthat::expect_s3_class(fm2,'species_mix')

})

