testthat::test_that('testing partial species mix bernoulli ', {

  library(ecomix)
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
  spp_form <- as.formula(~1+w1+w2)
  beta <- matrix(c(-3.6,0.5,
                   -0.9,1.0,
                   0.9,-2.9,
                   2.2,5.4),
                 4,2,byrow=TRUE)
  gamma <- matrix(c(rnorm(50,1),rnorm(50,-2)),50,2)
  dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form, spp_form, dat = dat,
                                         beta = beta, gamma = gamma,
                                         n_mixtures = 4,
                                         distribution = "bernoulli")

  archetype_formula <- sam_form
  species_formula <- spp_form

  test_dat <- ecomix:::clean_data_sam(simulated_data, archetype_formula, species_formula, distribution = 'negative_binomial')
  test_dat$mf.X
  test_dat$mf.W

  y <- simulated_data[,1:50]
  X <- ecomix:::get_X_sam(archetype_formula, test_dat$mf.X)
  W <- ecomix:::get_W_sam(species_formula, test_dat$mf.W)

  n_mixtures <- G <- 4
  S <- ncol(y)
  spp_weights <- rep(1,S)
  site_spp_weights <- matrix(1,nrow(y),S)
  disty <- 1
  y_is_na <- is.na(y)
  inits <- NULL
  control <- species_mix.control(quiet = FALSE)
  offset <- rep(0,nrow(X))

  fm_sp_mods <-  surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W,
                                       site_spp_weights, offset, y_is_na, disty,
                                       .parallel = control$cores, .verbose = FALSE)


  starting_values <- ecomix:::get_initial_values_sam(y = y, X = X, W = W,
                                                     spp_weights = spp_weights,
                                                     site_spp_weights = site_spp_weights,
                                                     offset = offset, y_is_na = y_is_na,
                                                     G = G, S = S,
                                                     disty = disty,
                                                     control = control)


  fits <- starting_values$fits
  taus <- starting_values$taus
  pis <- starting_values$pis
  first_fit <- starting_values$first_fit
  logls_mus <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)


  ss <- 1
  spp_conditional_max <- ecomix:::apply_glm_spp_coefs_sams(ss, y, X, W, G, taus,
                                                           site_spp_weights,
                                                           offset, y_is_na, disty, fits)

  fm_spp_coefs <- surveillance::plapply(seq_len(S),
                                        ecomix:::apply_glm_spp_coefs_sams,
                                        y, X, W, G, taus, site_spp_weights,
                                        offset, y_is_na, disty, fits,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

  gg <- 1
  mix_conditional_max <- ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W,
                                                           site_spp_weights,
                                                           offset, y_is_na, disty,
                                                           taus, fits, logls_mus$fitted)

  fm_mix_coefs <- surveillance::plapply(seq_len(G),
                                        ecomix:::apply_glm_mix_coefs_sams,
                                        y, X, W,
                                        site_spp_weights, offset, y_is_na, disty,
                                        taus, fits, logls_mus$fitted,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

  partial_ECM <- ecomix:::fitmix_ECM_sam(y, X, W, spp_weights, site_spp_weights,
                                         offset, y_is_na, G, S, disty,
                                         control=species_mix.control(em_steps=10))

  start_vals <- ecomix:::get_starting_values_sam(y = y, X = X, W = W,
                                                 spp_weights = spp_weights,
                                                 site_spp_weights = site_spp_weights,
                                                 offset = offset,
                                                 y_is_na = y_is_na,
                                                 G = G, S = S,
                                                 disty = disty,
                                                 control = species_mix.control(em_steps = 3))

  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control(optimise_cpp = 0,loglOnly_cpp = 1))


  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control(optimise_cpp = 0,derivOnly_cpp = 1))

  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control())

  test_part_sam <- species_mix(sam_form,spp_form,simulated_data,4,
                               distribution = 'bernoulli',
                               control = species_mix.control(em_steps = 5))

  object <- test_part_sam

  test_part_sam$vcov <- vcov(test_part_sam)

  summary(test_part_sam)

  species_mix.bootstrap(test_part_sam,nboot = 10)

})



testthat::test_that('testing partial species mix poisson ', {

  library(ecomix)
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
  spp_form <- as.formula(~1+w1+w2)
  alpha <- rnorm(50,-4,.5)
  beta <- matrix(c(-3.6,0.5,
                   -0.9,1.0,
                   0.9,-2.9,
                   2.2,5.4),
                 4,2,byrow=TRUE)
  gamma <- matrix(c(rnorm(50,1),rnorm(50,-2)),50,2)
  dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form, spp_form, dat = dat,
                                         alpha = alpha, beta = beta, gamma = gamma,
                                         n_mixtures = 4,
                                         distribution = "poisson")

  archetype_formula <- sam_form
  species_formula <- spp_form

  test_dat <- ecomix:::clean_data_sam(simulated_data, archetype_formula, species_formula, distribution = 'poisson')
  test_dat$mf.X
  test_dat$mf.W

  y <- simulated_data[,1:50]
  X <- ecomix:::get_X_sam(archetype_formula, test_dat$mf.X)
  W <- ecomix:::get_W_sam(species_formula, test_dat$mf.W)

  n_mixtures <- G <- 4
  S <- ncol(y)
  spp_weights <- rep(1,S)
  site_spp_weights <- matrix(1,nrow(y),S)
  disty <- 2
  y_is_na <- is.na(y)
  inits <- NULL
  control <- species_mix.control(quiet = FALSE)
  offset <- rep(0,nrow(X))

  fm_sp_mods <-  surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W,
                                       site_spp_weights, offset, y_is_na, disty,
                                       .parallel = control$cores, .verbose = FALSE)


  starting_values <- ecomix:::get_initial_values_sam(y = y, X = X, W = W,
                                                     spp_weights = spp_weights,
                                                     site_spp_weights = site_spp_weights,
                                                     offset = offset, y_is_na = y_is_na,
                                                     G = G, S = S,
                                                     disty = disty,
                                                     control = control)


  fits <- starting_values$fits
  taus <- starting_values$taus
  pis <- starting_values$pis
  first_fit <- starting_values$first_fit
  logls_mus <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)


  ss <- 1
  spp_conditional_max <- ecomix:::apply_glm_spp_coefs_sams(ss, y, X, W, G, taus,
                                                           site_spp_weights,
                                                           offset, y_is_na, disty, fits)

  fm_spp_coefs <- surveillance::plapply(seq_len(S),
                                        ecomix:::apply_glm_spp_coefs_sams,
                                        y, X, W, G, taus, site_spp_weights,
                                        offset, y_is_na, disty, fits,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

  gg <- 1
  mix_conditional_max <- ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W,
                                                           site_spp_weights,
                                                           offset, y_is_na, disty,
                                                           taus, fits, logls_mus$fitted)

  fm_mix_coefs <- surveillance::plapply(seq_len(G),
                                        ecomix:::apply_glm_mix_coefs_sams,
                                        y, X, W,
                                        site_spp_weights, offset, y_is_na, disty,
                                        taus, fits, logls_mus$fitted,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

  partial_ECM <- ecomix:::fitmix_ECM_sam(y, X, W, spp_weights, site_spp_weights,
                                         offset, y_is_na, G, S, disty,
                                         control=species_mix.control(em_steps=10))

  start_vals <- ecomix:::get_starting_values_sam(y = y, X = X, W = W,
                                                 spp_weights = spp_weights,
                                                 site_spp_weights = site_spp_weights,
                                                 offset = offset,
                                                 y_is_na = y_is_na,
                                                 G = G, S = S,
                                                 disty = disty,
                                                 control = species_mix.control(em_steps = 3))

  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control(optimise_cpp = 0,loglOnly_cpp = 1))


  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control(optimise_cpp = 0,derivOnly_cpp = 1))

  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control())

  test_part_sam <- species_mix(sam_form,spp_form,simulated_data,4,
                               distribution = 'poisson',
                               control = species_mix.control(em_steps = 5))
})

testthat::test_that('testing partial species mix ippm ', {

  library(ecomix)
  set.seed(123)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
  spp_form <- as.formula(~1+w1+w2)
  alpha <- rnorm(50,-13,.5)
  beta <- matrix(c(-2.6,-1.5,
                   -1.1,1.8,
                   0.9,-2.9,
                   2.2,1.4),
                 4,2,byrow=TRUE)
  gamma <- matrix(c(rnorm(50,-4,.5),rnorm(50,-2,.5)),50,2)
  dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
  dat[,-1] <- scale(dat[,-1])
  archetype_formula <- sam_form
  species_formula <- spp_form

  simulated_data <- species_mix.simulate(sam_form, spp_form, dat = dat,
                                         alpha = alpha, beta = beta, gamma = gamma,
                                         n_mixtures = 4,
                                         distribution = "ippm")

  test_dat <- ecomix:::clean_data_sam(simulated_data, archetype_formula, species_formula, distribution = 'ippm')
  test_dat$mf.X
  test_dat$mf.W

  y <- simulated_data[,1:50]
  X <- ecomix:::get_X_sam(archetype_formula, test_dat$mf.X)
  W <- ecomix:::get_W_sam(species_formula, test_dat$mf.W)

  n_mixtures <- G <- 4
  S <- ncol(y)
  spp_weights <- rep(1,S)
  site_spp_weights <- attr(simulated_data,'ippm_weights')
  disty <- 3
  y_is_na <- is.na(y)
  inits <- NULL
  control <- species_mix.control(quiet = FALSE)
  offset <- rep(0,nrow(X))

  fm_sp_mods <-  surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W,
                                       site_spp_weights, offset, y_is_na, disty,
                                       .parallel = control$cores, .verbose = FALSE)


  starting_values <- ecomix:::get_initial_values_sam(y = y, X = X, W = W,
                                                     spp_weights = spp_weights,
                                                     site_spp_weights = site_spp_weights,
                                                     offset = offset, y_is_na = y_is_na,
                                                     G = G, S = S,
                                                     disty = disty,
                                                     control = control)


  fits <- starting_values$fits
  taus <- starting_values$taus
  pis <- starting_values$pis
  first_fit <- starting_values$first_fit
  logls_mus <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)


  ss <- 1
  spp_conditional_max <- ecomix:::apply_glm_spp_coefs_sams(ss, y, X, W, G, taus,
                                                           site_spp_weights,
                                                           offset, y_is_na, disty, fits)

  fm_spp_coefs <- surveillance::plapply(seq_len(S),
                                        ecomix:::apply_glm_spp_coefs_sams,
                                        y, X, W, G, taus, site_spp_weights,
                                        offset, y_is_na, disty, fits,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

  gg <- 1
  mix_conditional_max <- ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W,
                                                           site_spp_weights,
                                                           offset, y_is_na, disty,
                                                           taus, fits, logls_mus$fitted)

  fm_mix_coefs <- surveillance::plapply(seq_len(G),
                                        ecomix:::apply_glm_mix_coefs_sams,
                                        y, X, W,
                                        site_spp_weights, offset, y_is_na, disty,
                                        taus, fits, logls_mus$fitted,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

  partial_ECM <- ecomix:::fitmix_ECM_sam(y, X, W, spp_weights, site_spp_weights,
                                         offset, y_is_na, G, S, disty,
                                         control=species_mix.control(em_steps=5))

  start_vals <- ecomix:::get_starting_values_sam(y = y, X = X, W = W,
                                                 spp_weights = spp_weights,
                                                 site_spp_weights = site_spp_weights,
                                                 offset = offset,
                                                 y_is_na = y_is_na,
                                                 G = G, S = S,
                                                 disty = disty,
                                                 control = species_mix.control(em_steps = 3))

  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control(optimise_cpp = 0,loglOnly_cpp = 1))


  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control(optimise_cpp = 0,derivOnly_cpp = 1))

  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,start_vals,
                               control=ecomix:::species_mix.control())

  test_part_sam <- species_mix(sam_form,spp_form,simulated_data,4,
                               weights = attr(simulated_data,'ippm_weights'),
                               distribution = 'ippm',
                               control = species_mix.control(em_steps = 5))


})

testthat::test_that('testing partial species mix negative binomial ', {

library(ecomix)
set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
spp_form <- as.formula(~1+w1+w2)
beta <- matrix(c(-3.6,0.5,
                 -0.9,1.0,
                  0.9,-2.9,
                  2.2,5.4),
                4,2,byrow=TRUE)
gamma <- matrix(c(rnorm(50,1),rnorm(50,-2)),50,2)
dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- species_mix.simulate(sam_form, spp_form, dat = dat,
                                       beta = beta, gamma = gamma,
                                       n_mixtures = 4,
                                       distribution = "negative_binomial")

archetype_formula <- sam_form
species_formula <- spp_form

test_dat <- ecomix:::clean_data_sam(simulated_data, archetype_formula, species_formula, distribution = 'negative_binomial')
test_dat$mf.X
test_dat$mf.W

y <- simulated_data[,1:50]
X <- ecomix:::get_X_sam(archetype_formula, test_dat$mf.X)
W <- ecomix:::get_W_sam(species_formula, test_dat$mf.W)
n_mixtures <- G <- 4
S <- ncol(y)
spp_weights <- rep(1,S)
site_spp_weights <- matrix(1,nrow(y),S)
disty <- 4
y_is_na <- is.na(y)
inits <- NULL
control <- species_mix.control(quiet = FALSE)
offset <- rep(0,nrow(X))
fm_sp_mods <-  surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W,
                                     site_spp_weights, offset, y_is_na, disty,
                                     .parallel = control$cores, .verbose = FALSE)


starting_values <- ecomix:::get_initial_values_sam(y = y, X = X, W = W,
                                                   spp_weights = spp_weights,
                                                   site_spp_weights = site_spp_weights,
                                                   offset = offset, y_is_na = y_is_na,
                                                   G = G, S = S,
                                                   disty = disty,
                                                   control = control)

fits <- starting_values$fits
taus <- starting_values$taus
pis <- starting_values$pis
first_fit <- starting_values$first_fit
logls_mus <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)


ss <- 1
spp_conditional_max <- ecomix:::apply_glm_spp_coefs_sams(ss, y, X, W, G, taus,
                                                         site_spp_weights,
                                                         offset, y_is_na, disty, fits)

fm_spp_coefs <- surveillance::plapply(seq_len(S),
                                      ecomix:::apply_glm_spp_coefs_sams,
                                      y, X, W, G, taus, site_spp_weights,
                                      offset, y_is_na, disty, fits,
                                      .parallel = control$cores,
                                      .verbose = FALSE)

tmp <- nlminb(start=fits$beta, objective=ecomix:::incomplete_negbin_logl, gradient=NULL,
              hessian=NULL, pis=pis, first_fit=first_fit, fits=fits, G=G, S=S)

partial_ECM <- ecomix:::fitmix_ECM_sam(y, X, W, spp_weights, site_spp_weights,
                                       offset, y_is_na, G, S, disty,
                                       control=species_mix.control(em_steps=10))

start_vals <- ecomix:::get_starting_values_sam(y = y, X = X, W = W,
                                      spp_weights = spp_weights,
                                      site_spp_weights = site_spp_weights,
                                      offset = offset,
                                      y_is_na = y_is_na,
                                      G = G, S = S,
                                      disty = disty,
                                      control = species_mix.control(em_steps = 3))

tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                             S,G,disty,start_vals,
                             control=ecomix:::species_mix.control(optimise_cpp = 0,loglOnly_cpp = 1))


tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                             S,G,disty,start_vals,
                             control=ecomix:::species_mix.control(optimise_cpp = 0,derivOnly_cpp = 1))

start_vals$beta[]<-1

tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                             S,G,disty,start_vals,
                             control=ecomix:::species_mix.control())

test_part_sam <- species_mix(sam_form,spp_form,simulated_data,4,
                             distribution = 'negative_binomial',
                             control = species_mix.control(em_steps = 5))

})

testthat::test_that('testing partial species mix S3 classes', {

  AIC()

  BIC()

  coef()

  preds <- predict()


})

