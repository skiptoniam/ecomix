context('species_mix ippm')

library(ecomix)
library(raster)
library(scales)

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
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,thetas,dist="ippm")
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
  one_sp_ippm <- ecomix:::apply_glm_sam_inits(ss = ss, y = y, X = X, site_spp_weights = site_spp_weights, offset = offset, y_is_na = y_is_na,disty = disty)
  testthat::expect_is(one_sp_ippm,'list')

  # check that many species ippms work - expect back a list.
  all_sp_ippm <-surveillance::plapply(seq_len(S), ecomix:::apply_glm_sam_inits, y, X, site_spp_weights, offset, y_is_na,disty)
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

  # tmp <- ecomix:::get_starting_values_sam(y,X,spp_weights,site_spp_weights,offset,y_is_na,G,S,disty,control)

  # does a ippm work?
  # wts <- rbind(presence_sites[,-1],background_sites[,-1])
  # colnames(wts) <- c(sp_name)#,"const","x1","x2")
  # wts <- as.matrix(wts)
  offset <- rep(0,nrow(as.matrix(site_spp_weights)))
  sam_form <- as.formula(paste0('cbind(',paste(colnames(y),collapse = ','),")~1+x1+x2"))
  sp_form <- ~ 1
  model_data <- make_mixture_data(y,X[,-1])

  # print(head(site_spp_weights))
  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'ippm',
                     weights = as.matrix(site_spp_weights),
                     n_mixtures = 4,
                     control = species_mix.control(minimum_sites_prevelance = 50,init_method = 'kmed'))
  # testthat::expect_s3_class(fm1,'ippm')
  testthat::expect_s3_class(fm1,'species_mix')

  # expect error if there are no weights
  testthat::expect_error(  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'ippm',weights = 1,
                                              n_mixtures=4))

  fm3 <- species_mix(sam_form, sp_form, model_data, distribution = 'ippm',weights = as.matrix(site_spp_weights), n_mixtures=4)
  # test if species weight have different names.
  colnames(site_spp_weights) <- NULL#paste0('blah',seq_len(S))
  testthat::expect_error(fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'ippm',
                     weights = as.matrix(site_spp_weights),
                     n_mixtures=4))
})
