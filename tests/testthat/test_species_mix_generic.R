context('species_mix generic functions')
library(ecomix)

testthat::test_that('species mix internal functions classes work', {


  set.seed(42)
  archetype_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~x1+x2"))
  theta <- matrix(c(1,-2.9,-3.6,1,-0.9,1,1,.9,7.9),3,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(archetype_formula = archetype_form, species_formula = ~1,dat = dat,dist="bernoulli")
  model_data <- simulated_data

  #test formula error
  testthat::expect_error(fm1 <- species_mix(NA, ~1, model_data, distribution = 'bernoulli', n_mixtures=3))
  testthat::expect_error(fm3 <- species_mix(NA, ~1, model_data, distribution = 'poisson', n_mixtures=3))
  testthat::expect_error(fm4 <- species_mix(NA, ~1, model_data, distribution = 'ippm', n_mixtures=3))
  testthat::expect_error(fm5 <- species_mix(NA, ~1, model_data, distribution = 'negative_binomial', n_mixtures=3))

  ## test the internal functions
  ## test to see if the species formula checks are working.

  f1 <- y ~ x
  f2 <- y ~ x + z
  f3 <- y ~ 1
  f4 <- NULL

  testthat::expect_true(ecomix:::check_species_formula(f1)==2)
  testthat::expect_true(ecomix:::check_species_formula(f2)==2)
  testthat::expect_true(ecomix:::check_species_formula(f3)==1)
  testthat::expect_true(ecomix:::check_species_formula(f4)==0)

  y <- simulated_data[,1:20]
  X <- simulated_data[,-1:-21]
  W <- simulated_data[,21,drop=FALSE]
  offset <- rep(0,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- length(attr(simulated_data,"pis"))
  S <- 20
  nPX <- ncol(X)
  nPW <- ncol(W)
  control <- species_mix.control()

  # test a new glm function bernoulli
  ss <- 1
  disty <- 1
  fm1 <- ecomix:::apply_glmnet_sam_inits(ss, y, X, W, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,4)

  fm_bern <- surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)

  alphas <- lapply(fm_bern, `[[`, 1)
  testthat::expect_length(unlist(alphas),S)

  betas <- lapply(fm_bern, `[[`, 2)
  testthat::expect_length(do.call(rbind, betas),S*nPX)

  gammas <- lapply(fm_bern, `[[`, 3)
  testthat::expect_length(do.call(rbind, gammas),S)

  thetas <- unlist(lapply(fm_bern, `[[`, 4))
  testthat::expect_true(all(thetas==-99999))

  ## poisson
  simulated_data <- species_mix.simulate(sam_form,~1,dat,dist="poisson")

  y <- simulated_data[,1:20]
  ss <- 1
  disty <- 2
  fm1 <- ecomix:::apply_glmnet_sam_inits(ss, y, X, W, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,4)

  fm_pois <- surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)

  alphas <- lapply(fm_pois, `[[`, 1)
  testthat::expect_length(unlist(alphas),S)

  betas <- lapply(fm_pois, `[[`, 2)
  testthat::expect_length(do.call(rbind, betas),S*nPX)

  gammas <- lapply(fm_pois, `[[`, 3)
  testthat::expect_length(do.call(rbind, gammas),S)

  thetas <- unlist(lapply(fm_pois, `[[`, 4))
  testthat::expect_true(all(thetas==-99999))
  })


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
                                         dat,beta=beta,dist="bernoulli")
  testthat::expect_message(fm1 <- species_mix(NULL, sp_form,
                                              model_data,
                                              distribution = 'bernoulli',
                                              n_mixtures=3))

  dup_spp_data <- cbind('spp1'=model_data[,1],model_data)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',c(1,1:20)),collapse = ','),")~1+x1+x2"))

  testthat::expect_message(fm1 <- species_mix(sam_form, sp_form, dup_spp_data,
                                              distribution = 'bernoulli',
                                              n_mixtures=3))

})

testthat::test_that('testing species mix S3 class functions', {

  library(ecomix)
  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',
                                       paste(paste0('spp',1:20),collapse = ','),
                                       ")~1+x1+x2"))
  sp_form <- ~ 1
  beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,1.9),3,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  model_data <- species_mix.simulate(archetype_formula=sam_form,
                                         species_formula=sp_form,
                                         dat, beta= beta,
                                         dist="bernoulli")
  fm1 <- species_mix(sam_form, sp_form, model_data,
                     distribution = 'bernoulli',
                     n_mixtures=3)

  print(fm1)
  testthat::expect_error(summary(fm1))
  fm1$vcov <- vcov(fm1)
  testthat::expect_is(summary(fm1),'matrix')
  testthat::expect_length(AIC(fm1),1)
  testthat::expect_length(BIC(fm1),1)
  testthat::expect_is(coef(fm1),'list')

})


testthat::test_that('species mix generic vcov functions', {

  # build and test a single model.
  # estimate variance-covariance matrix
  library(ecomix)
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:100),collapse = ','),")~1+x1+x2"))
  theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm1 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'bernoulli', n_mixtures=4)

  vcv_mat <- vcov(object = fm1)
  testthat::expect_equal(nrow(vcv_mat),nrow(vcv_mat))
  testthat::expect_is(vcv_mat,'matrix')
  testthat::expect_true(all(is.finite(sqrt(diag(vcv_mat)))))

  vcv_mat_bb <- vcov(object = fm1,method = 'BayesBoot', nboot = 10)
  testthat::expect_equal(nrow(vcv_mat),nrow(vcv_mat))
  testthat::expect_is(vcv_mat_bb,'matrix')
  testthat::expect_true(all(is.finite(sqrt(diag(vcv_mat_bb)))))

})

testthat::test_that('species mix predict functions', {

  # build and test a single model.
  library(ecomix)
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:100),collapse = ','),")~1+x1+x2"))
  theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm1 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'bernoulli', n_mixtures=4)

  preds <- predict(fm1)
  testthat::expect_length(preds,2)
  testthat::expect_is(preds,'list')

  dat2 <- data.frame(x1=runif(100,-2.5,2.5),x2=rnorm(100,-2.5,2.5))
  preds2 <- predict(fm1, newobs = dat2)
  testthat::expect_is(preds2,'list')

  # poisson
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="poisson")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm2 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'poisson', n_mixtures=4)

  preds3 <- predict(fm2)
  testthat::expect_length(preds3,2)
  testthat::expect_is(preds3,'list')

  preds4 <- predict(fm2, newobs = dat2)
  testthat::expect_is(preds4,'list')

  # negative binomial
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="negative_binomial")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm3 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'negative_binomial', n_mixtures=4)

  preds5 <- predict(fm3)
  testthat::expect_length(preds3,2)
  testthat::expect_is(preds5,'list')

  preds6 <- predict(fm2, newobs = dat2)
  testthat::expect_is(preds6,'list')

  #gaussian
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="gaussian")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm4 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'gaussian', n_mixtures=4)

  preds7 <- predict(fm4)
  testthat::expect_length(preds7,2)
  testthat::expect_is(preds7,'list')

  preds8 <- predict(fm4, newobs = dat2)
  testthat::expect_is(preds8,'list')

  testthat::expect_error(preds8 <- predict('a'))

})


testthat::test_that('species mix gaussian', {

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  theta <- matrix(c(-2.9,3.6,0.5,1,-0.9,1,.9,4.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="gaussian")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])

  y <- simulated_data$species_data
  X <- simulated_data$covariate_data
  offset <- rep(0,nrow(y))
  weights <- rep(1,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- length(simulated_data$pi)
  S <- length(simulated_data$sp.int)
  control <- species_mix.control()
  disty <- 6

  # test a single gaussian model
  i <- 1
  # for(i in 1:S)
  testthat::expect_length(ecomix:::apply_glmnet_sam_inits(i, y, X, site_spp_weights, offset, y_is_na, disty),3)
  fm_gaussianint <- surveillance::plapply(1:S, ecomix:::apply_glmnet_sam_inits,  y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_gaussianint)[1,],S)

  #get the taus
  starting_values <- ecomix:::initiate_fit_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
  fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,disp=starting_values$disp)
  first_fit <- list(x = X, y = y, weights=site_spp_weights, offset=offset)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  taus <- ecomix:::shrink_taus(taus, max_tau=1/G + 0.1, G)

  ## get to this in a bit
  gg <- 1
  testthat::expect_length(ecomix:::apply_glm_group_tau_sam(gg, y, X, site_spp_weights, offset, y_is_na, disty, taus, fits, logls$fitted),2)

  # ## now let's try and fit the optimisation
  sv <- ecomix:::get_starting_values_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
  tmp <- ecomix:::sam_optimise(y,X,offset,spp_weights,site_spp_weights, y_is_na, S, G, disty, start_vals = sv, control)
  testthat::expect_length(tmp,17)

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  sp_form <- ~1
  theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="gaussian")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'gaussian',
                     n_mixtures=3)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'gaussian',
                     n_mixtures=3,control=species_mix.control(em_prefit = FALSE))
  testthat::expect_s3_class(fm2,'species_mix')

})


testthat::test_that('species mix poisson', {

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="poisson")
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
  disty <- 2

  # test a single poisson model
  i <- 1
  testthat::expect_length(ecomix:::apply_glmnet_sam_inits(i, y, X, site_spp_wts, offset, y_is_na, disty),3)
  fm_poissonint <- surveillance::plapply(1:S, ecomix:::apply_glmnet_sam_inits,  y, X, site_spp_wts, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_poissonint)[1,],S)

  # test that the starting values work.
  testthat::expect_length(tmp <- ecomix:::get_starting_values_sam(y, X, spp_wts, site_spp_wts, offset, y_is_na, G, S, disty, control),9)

  #get the taus
  starting_values <- ecomix:::initiate_fit_sam(y, X,spp_wts, site_spp_wts, offset, y_is_na, G, S, disty, control)
  fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,disp=starting_values$disp)
  first_fit <- list(x = X, y = y, weights=site_spp_wts, offset=offset)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(first_fit, fits, spp_wts, G, S, disty)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  taus <- ecomix:::shrink_taus(taus, max_tau=1/G + 0.1, G)

  ## get to this in a bit
  gg <- 1
  testthat::expect_length(ecomix:::apply_glm_group_tau_sam(gg, y, X, site_spp_wts, offset, y_is_na, disty, taus,fits,logls$fitted),2)

  # ## now let's try and fit the optimisation
  sv <- ecomix:::get_starting_values_sam(y, X, spp_wts, site_spp_wts, offset, y_is_na, G, S, disty, control)
  tmp <- ecomix:::sam_optimise(y,X,offset,spp_wts,site_spp_wts, y_is_na, S, G, disty, start_vals = sv, control)
  testthat::expect_length(tmp,17)

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  sp_form <- ~1
  theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="poisson")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'poisson',
                     n_mixtures=3)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'poisson',
                     n_mixtures=3,control=species_mix.control(em_prefit = FALSE))
  testthat::expect_s3_class(fm2,'species_mix')

})


context('species_mix bernoulli')
library(ecomix)


testthat::test_that('species mix bernoulli', {


  rm(list = ls())
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])

  y <- simulated_data$species_data
  X <- simulated_data$covariate_data[,-1]
  W <- simulated_data$covariate_data[,1,drop=FALSE]
  offset <- rep(0,nrow(y))
  weights <- rep(1,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- length(simulated_data$pi)
  S <- length(simulated_data$sp.int)
  control <- species_mix.control()
  disty <- 1

  # test a single bernoulli model
  i <- 1
  testthat::expect_length(ecomix:::apply_glmnet_sam_inits(i, y, X, W,
                                                       site_spp_weights, offset,
                                                       y_is_na, disty),4)
  fm_bernoulliint <- surveillance::plapply(1:S, ecomix:::apply_glmnet_sam_inits,
                                           y, X, W, site_spp_weights, offset, y_is_na,
                                           disty,
                                           .parallel = control$cores,
                                           .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_bernoulliint)[1,],S)

  #get the taus
  starting_values <- ecomix:::initiate_fit_sam(y, X, W, site_weights,
                                               site_spp_weights, offset,
                                               y_is_na, G, S, disty, control)
  fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,
               gamma=starting_values$gamma,
               theta=starting_values$theta)
  first_fit <- list(x = X, W = W, y = y, site_spp_weights=site_spp_weights,
                    offset=offset)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  taus <- ecomix:::shrink_taus(taus, max_tau=1/G + 0.1, G)

  ## get to this in a bit
  ss <- 1
  test <- ecomix:::apply_glmnet_spp_coefs_sams(ss, y, X, W, G, taus,
                                               site_spp_weights,
                                               offset, y_is_na,
                                               disty, fits)
  gg <- 1
  testthat::expect_length(ecomix:::apply_glm_group_tau_sam(gg, y, X,
                                                           site_spp_weights,
                                                           offset, y_is_na,
                                                           disty, taus, fits,
                                                           logls$fitted),2)

  # ## now let's try and fit the optimisation
  sv <- ecomix:::get_starting_values_sam(y, X, spp_weights, site_spp_weights,
                                         offset, y_is_na, G, S, disty, control)
  tmp <- ecomix:::sam_optimise(y,X,offset,spp_weights,site_spp_weights, y_is_na,
                               S, G, disty, start_vals = sv, control)
  testthat::expect_length(tmp,17)


  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))
  sp_form <- ~ 1
  theta <- matrix(c(1,-2.9,-3.6,
                    1,-0.9,1,
                    1,.9,7.9),
                  3,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(archetype_formula=sam_form, species_formula=sp_form,
                                              dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli',
                     n_mixtures=3)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli',
                     n_mixtures=3,control=species_mix.control(em_prefit = FALSE))
  testthat::expect_s3_class(fm2,'species_mix')


})


testthat::test_that('species mix ippm', {

  # library(ecomix)
  # library(raster)
  set.seed(42)

  x1 <- sort(runif(1000,-2.5,2.5))
  x2 <- I(x1)^2

  n_g <- 4
  n_sp <- 50
  set.seed(123)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  thetas <- matrix(c( 1, 4.0,-3.0,
                      1, 1.8,   0,
                      1,-1.2, 0.1,
                      1,-4.2,-4.8),4,3,byrow=TRUE)

  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,thetas,dist="ippm")
  # wts <- simulated_data$background_weights
  offset <- simulated_data$offset

  ## test the internal functions.
  ## test the apply functions for ippm

  y <- simulated_data$species_data
  X <- simulated_data$covariate_data
  y_is_na <- is.na(y)
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- as.matrix(simulated_data$background_weights)
  G <- 4
  S <- n_sp
  ss <- 1
  disty <- 3
  nP <- dim(thetas)[2]-1
  control <- species_mix.control(minimum_sites_occurrence = 50)

  # test if one species ippm working - expect matrix of coefs back

  one_sp_ippm <- ecomix:::apply_glmnet_sam_inits(ss = ss, y = y, X = X, site_spp_weights = site_spp_weights, offset = offset, y_is_na = y_is_na,disty = disty)
  testthat::expect_is(one_sp_ippm,'list')

  # check that many species ippms work - expect back a list.
  all_sp_ippm <-surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, site_spp_weights, offset, y_is_na,disty)
  testthat::expect_is(all_sp_ippm,'list')

  alpha <- lapply(all_sp_ippm, `[[`, 1)
  testthat::expect_length(unlist(alpha),S)

  beta <- lapply(all_sp_ippm, `[[`, 2)
  testthat::expect_length(do.call(rbind, beta),S*nP)

  disp <- unlist(lapply(all_sp_ippm, `[[`, 3))
  testthat::expect_true(all(is.na(disp)))

  #glmnet coefs
  all_coefs_mat <- t(do.call(cbind,beta))
  # mix_coefs <- all_coefs_mat[,-1] # drop intercepts
  tmp1 <- kmeans(all_coefs_mat, centers=G, nstart=100)
  tmp_grp <- tmp1$cluster
  grp_coefs <- apply(all_coefs_mat, 2, function(x) tapply(x, tmp_grp, mean))

  #now we need to estimate the taus.
  S <- 50
  G <- 4

  # expect error if wrong data is in the starting values
  testthat::expect_error(  starting_values <- ecomix:::initiate_fit_sam(NULL, X, weights, offset, y_is_na, G, S, control))
  testthat::expect_error(  starting_values <- ecomix:::initiate_fit_sam(y, NULL, weights, offset, y_is_na, G, S, control))
  testthat::expect_error(  starting_values <- ecomix:::initiate_fit_sam(y, X, NULL, offset, y_is_na, G, S, control))
  testthat::expect_error(  starting_values <- ecomix:::initiate_fit_sam(y, X, weights, offset, NULL, G, S, control))

  #expect list back
  starting_values <- ecomix:::initiate_fit_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
  testthat::expect_is(starting_values,'list')

  fits <- list(beta=starting_values$beta, alpha=starting_values$alpha)
  first_fit <- list(x = X, y = y, site_spp_weights=site_spp_weights, offset=offset, y_is_na=y_is_na)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  testthat::expect_is(logls$logl_sp,'matrix')
  testthat::expect_equal(ncol(logls$logl_sp), G)
  testthat::expect_equal(nrow(logls$logl_sp), S)

  # estimate the posteriors for taus
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  testthat::expect_is(taus,'matrix')
  testthat::expect_equal(ncol(taus), G)
  testthat::expect_equal(nrow(taus), S)

  # skrink the taus
  taus <- ecomix:::shrink_taus(taus, max_tau=0.8, G)

  ## now test if the group_tau glm works
  testthat::expect_error(fm_g1 <- ecomix:::apply_glm_group_tau_sam(gg = 1, y = y, X = X, weights = weights, offset = offset,
                                                                   y_is_na = y_is_na, tau = tau, return_all_coefs = FALSE))
  testthat::expect_error(fm_g1 <- ecomix:::apply_glm_group_tau_sam(gg = 1, y = NULL, X = X, weights = weights, offset = offset,
                                                                   y_is_na = y_is_na, tau = taus, return_all_coefs = FALSE))
  testthat::expect_error(fm_g1 <- ecomix:::apply_glm_group_tau_sam(gg = 1, y = y, X = NULL, weights = weights, offset = 'blah',
                                                                   y_is_na = y_is_na, tau = taus, return_all_coefs = FALSE))
  testthat::expect_error(fm_g1 <- ecomix:::apply_glm_group_tau_sam(gg = 1, y = y, X = X, weights = 'a', offset = offset,
                                                                   y_is_na = y_is_na, tau = taus, return_all_coefs = FALSE))
  fm_g1 <- ecomix:::apply_glm_group_tau_sam(gg = 1, y = y, X = X, site_spp_weights = site_spp_weights, offset=offset, y_is_na = y_is_na,
                                            disty = disty,  tau = taus, fits, logls_mus$fitted)

  testthat::expect_is(fm_g1,'numeric')

  all_grp_ippm1 <- surveillance::plapply(seq_len(G), ecomix:::apply_glm_group_tau_sam, y, X, site_spp_weights, offset, y_is_na, disty, taus, fits, logls_mus$fitted)
  beta <- t(do.call(cbind,all_grp_ippm1))
  testthat::expect_is(beta,'matrix')

  # does a ippm work?
  offset <- rep(0,nrow(as.matrix(site_spp_weights)))
  sam_form <- as.formula(paste0('cbind(',paste(colnames(y),collapse = ','),")~1+x1+x2"))
  sp_form <- ~ 1
  model_data <- make_mixture_data(y,X[,-1])

  # print(head(site_spp_weights))
  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'ippm',
                     weights = as.matrix(site_spp_weights),
                     n_mixtures = 4,
                     control = species_mix.control(em_prefit = FALSE, minimum_sites_prevelance = 50,init_method = 'kmed'))
  testthat::expect_s3_class(fm1,'species_mix')

  # expect error if there are no weights
  testthat::expect_error(  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'ippm',weights = 1,
                                              n_mixtures=4))

  # test if species weight have different names.
  colnames(site_spp_weights) <- NULL#paste0('blah',seq_len(S))
  testthat::expect_error(fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'ippm',
                                            weights = as.matrix(site_spp_weights),
                                            n_mixtures=4))
})


context('species_mix negative binomial')
library(ecomix)



testthat::test_that('species_mix negative binomial', {

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2+z1+z2"))
  sp_form <- NULL
  theta <- matrix(c(1,-2.1, 1.6,  0.5,-.1,
                    1,-0.9,  -1, -4.3, 2,
                    1, 1.9, 3.9,  0.3,-2),3,5,byrow=TRUE)
  theta <- matrix(c(1,1.6,0.5, 0.5,-.1,
                    1,-7.9,3.6, -4.3, 2,
                    1,4.9,-2.9, 0.3,-2,
                    1,-0.2,-0.4,1,-3),4,5,byrow=TRUE)
  x1<-runif(100,0,2.5)
  z1<-rnorm(100,0,2.5)
  dat <- data.frame(y=rep(1,100),x1,x2=x1^2,z1,z2=z1^2)
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(archetype_formula=sam_form,
                                              species_formula=sp_form,
                                              dat,theta,dist="negative_binomial")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])

  y <- simulated_data$species_data
  X <- simulated_data$covariate_data
  offset <- rep(0,nrow(y))
  weights <- rep(1,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- length(simulated_data$pi)
  S <- length(simulated_data$sp.int)
  control <- species_mix.control()
  disty <- 4
  nP <- 4

  # test_nb <- nbglm::glm.fit.nbinom(X,y[,3],offset,weights,est_var = FALSE)
  # test_nb <- nbglm::glm.fit.nbinom(X,y[,1],offset,weights,est_var = TRUE)

  ss <- 1
  fm1 <- ecomix:::apply_glmnet_sam_inits(ss, y, X, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,3)
  #
  fm_nb <- surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  #
  alpha <- lapply(fm_nb, `[[`, 1)
  testthat::expect_length(unlist(alpha),S)

  beta <- lapply(fm_nb, `[[`, 2)
  testthat::expect_length(do.call(rbind, beta),S*nP)

  disp <- unlist(lapply(fm_nb, `[[`, 3))
  testthat::expect_length(disp,S)


  # # test a single negative_binomial model
  i <- 1
  testthat::expect_length(ecomix:::apply_glmnet_sam_inits(i, y, X, site_spp_weights, offset, y_is_na, disty),3)
  fm_negative_binomialint <- surveillance::plapply(1:S, ecomix:::apply_glmnet_sam_inits,  y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_negative_binomialint)[1,],S)

  # test that the starting values work.
  # testthat::expect_length(tmp <- ecomix:::fitmix_ECM_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control),8)

  #get the taus
  starting_values <- ecomix:::initiate_fit_sam(y, X, spp_weights,
                                               site_spp_weights, offset,
                                               y_is_na, G, S, disty, control)
  fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,disp=starting_values$disp)
  first_fit <- list(x = X, y = y, weights=site_spp_weights, offset=offset)

  # # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls$logl_sp, G, S)
  taus <- ecomix:::shrink_taus(taus, max_tau=1/G + 0.1, G)

  # ## get to this in a bit
  gg <- 1
  testthat::expect_length(ecomix:::apply_glm_group_tau_sam(gg, y, X, site_spp_weights, offset, y_is_na, disty, taus, fits, logls$fitted),4)
  #
  # # ## now let's try and fit the optimisation
  # start_vals <- ecomix:::get_starting_values_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
  #
  # tmp <- ecomix:::sam_optimise(y,X,offset,spp_weights,site_spp_weights, y_is_na, S, G, disty, start_vals, control)
  # testthat::expect_length(tmp,17)

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  sp_form <- ~1
  theta <- matrix(c(1,1.6,0.5,
                    1,-7.9,3.6,
                    1,4.9,-2.9,
                    1,-0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="negative_binomial")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  # fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'negative_binomial',
  #                    n_mixtures=4,control=species_mix.control(em_refit = 5, em_steps = 10))
  # testthat::expect_s3_class(fm1,'species_mix')
  # #
  # fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'negative_binomial',
  #                    n_mixtures=4,control=species_mix.control(em_prefit = FALSE))
  # testthat::expect_s3_class(fm2,'species_mix')
})


testthat::test_that('species mix one covariate and one group', {

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  theta <- matrix(c(1,-1.6,-2.5,
                    1,-0.9,1,
                    1,2.9,2.9,
                    1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="bernoulli")
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
  simulated_data <- species_mix.simulate(sam_form,~1,dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])

  sp_form <- ~ 1
  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli',
                     # weights = as.matrix(site_spp_weights),
                     n_mixtures = 1,
                     control = species_mix.control(minimum_sites_prevelance = 50,init_method = 'kmed'))

})
