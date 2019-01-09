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
  simulated_data <- simulate_species_mix_data(archetype_formula=sam_form,
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
  fm1 <- ecomix:::apply_glm_sam_inits(ss, y, X, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,3)
  #
  fm_nb <- surveillance::plapply(seq_len(S), ecomix:::apply_glm_sam_inits, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  #
  alpha <- lapply(fm_nb, `[[`, 1)
  testthat::expect_length(unlist(alpha),S)

  beta <- lapply(fm_nb, `[[`, 2)
  testthat::expect_length(do.call(rbind, beta),S*nP)

  disp <- unlist(lapply(fm_nb, `[[`, 3))
  testthat::expect_length(disp,S)


  # # test a single negative_binomial model
  i <- 1
  testthat::expect_length(ecomix:::apply_glm_sam_inits(i, y, X, site_spp_weights, offset, y_is_na, disty),3)
  fm_negative_binomialint <- surveillance::plapply(1:S, ecomix:::apply_glm_sam_inits,  y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_negative_binomialint)[1,],S)

  # test that the starting values work.
  testthat::expect_length(tmp <- ecomix:::fitmix_EM_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control),8)

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
  start_vals <- ecomix:::get_starting_values_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
  #
  tmp <- ecomix:::sam_optimise(y,X,offset,spp_weights,site_spp_weights, y_is_na, S, G, disty, start_vals, control)
  testthat::expect_length(tmp,17)

  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  sp_form <- ~1
  theta <- matrix(c(1,1.6,0.5,
                    1,-7.9,3.6,
                    1,4.9,-2.9,
                    1,-0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="negative_binomial")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'negative_binomial',
                      n_mixtures=4,control=species_mix.control(em_refit = 5, em_steps = 10))
  testthat::expect_s3_class(fm1,'species_mix')
  #
  fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'negative_binomial',
                     n_mixtures=4,control=species_mix.control(em_prefit = FALSE))
  testthat::expect_s3_class(fm2,'species_mix')
})
