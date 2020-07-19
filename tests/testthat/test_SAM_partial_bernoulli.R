testthat::test_that('testing partial species mix bernoulli ', {

  library(ecomix)
  set.seed(42)
  nspp <- 20
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:nspp),
                                                      collapse = ','),
                                       ")~1+x1+x2"))
  sp_form <- as.formula(~1+w1+w2)
  beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,7.9),3,2,byrow=TRUE)
  gamma <- matrix(c(rnorm(50,1),rnorm(50,-2)),nspp,2)
  dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(archetype_formula=sam_form,
                                     species_formula=sp_form,
                                     dat,beta=beta,gamma = gamma,
                                     dist="bernoulli")

  archetype_formula <- sam_form
  species_formula <- sp_form

  test_dat <- ecomix:::clean_data_sam(simulated_data, archetype_formula, species_formula, distribution = 'negative_binomial')
  test_dat$mf.X
  test_dat$mf.W

  y <- simulated_data[,1:nspp]
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
  size <- rep(1,nrow(y))

  fm_sp_mods <-  surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W,
                                       site_spp_weights, offset, y_is_na, disty, size,
                                       .parallel = control$cores, .verbose = FALSE)


  starting_values <- ecomix:::get_initial_values_sam(y = y, X = X, W = W,
                                                     spp_weights = spp_weights,
                                                     site_spp_weights = site_spp_weights,
                                                     offset = offset, y_is_na = y_is_na,
                                                     G = G, S = S,
                                                     disty = disty,size = size,
                                                     control = control)


  fits <- starting_values$fits
  taus <- starting_values$taus
  pis <- starting_values$pis
  first_fit <- starting_values$first_fit
  logls_mus <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)


  ss <- 1
  spp_conditional_max <- ecomix:::apply_glm_spp_coefs_sams(ss, y, X, W, G, taus,
                                                           site_spp_weights,
                                                           offset, y_is_na, disty, fits, size)

  fm_spp_coefs <- surveillance::plapply(seq_len(S),
                                        ecomix:::apply_glm_spp_coefs_sams,
                                        y, X, W, G, taus, site_spp_weights,
                                        offset, y_is_na, disty, fits, size,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

  gg <- 1
  mix_conditional_max <- ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W,
                                                           site_spp_weights,
                                                           offset, y_is_na, disty,
                                                           taus, fits, logls_mus$fitted, size)

  fm_mix_coefs <- surveillance::plapply(seq_len(G),
                                        ecomix:::apply_glm_mix_coefs_sams,
                                        y, X, W,
                                        site_spp_weights, offset, y_is_na, disty,
                                        taus, fits, logls_mus$fitted, size,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

  partial_ECM <- ecomix:::fitmix_ECM_sam(y, X, W, spp_weights, site_spp_weights,
                                         offset, y_is_na, G, S, disty, size,
                                         control=species_mix.control(em_steps=10))

  start_vals <- ecomix:::get_starting_values_sam(y = y, X = X, W = W,
                                                 spp_weights = spp_weights,
                                                 site_spp_weights = site_spp_weights,
                                                 offset = offset,
                                                 y_is_na = y_is_na,
                                                 G = G, S = S,
                                                 disty = disty, size = size,
                                                 control = species_mix.control(em_steps = 10))

  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,size,start_vals,
                               control=ecomix:::species_mix.control(optimise_cpp = 0,loglOnly_cpp = 1))
  tmp[[1]]


  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,size,start_vals,
                               control=ecomix:::species_mix.control(optimise_cpp = 0,derivOnly_cpp = 1))
  tmp$scores
  tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                               S,G,disty,size,start_vals,
                               control=ecomix:::species_mix.control())

  test_part_sam <- species_mix(sam_form,sp_form,simulated_data,4,
                               distribution = 'bernoulli',
                               control = species_mix.control(em_steps = 5))

  test_part_sam <- species_mix(sam_form,sp_form,simulated_data,4,
                               distribution = 'bernoulli',
                               control = species_mix.control(em_steps = 2,
                                                             getscores_cpp = TRUE))
  fm <- species_mix(sam_form,sp_form,simulated_data,4,
                               distribution = 'bernoulli',
                               standardise = TRUE,
                               # inits = unlist(coef(test_part_sam)),
                               control = species_mix.control(em_steps = 2,
                                                             getscores_cpp = TRUE))

  fm <- species_mix(sam_form,sp_form,simulated_data,4,
                    distribution = 'bernoulli',
                    # standardise = TRUE,
                    inits = unlist(test_part_sam$coefs),
                    control = species_mix.control(em_steps = 2,
                                                  getscores_cpp = TRUE))


  test_part_sam$vcov <- vcov(test_part_sam,method = "BayesBoot", nboot=10)
  summary(test_part_sam)

  boots <- species_mix.bootstrap(test_part_sam, nboot = 10)
  preds <- predict(test_part_sam,boots)
  testthat::expect_is(preds,'list')

})
