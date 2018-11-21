context('species_mix negative binomial')
library(ecomix)



testthat::test_that('species_mix negative binomial', {


  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2+z1+z2"))
  sp_form <- NULL
  theta <- matrix(c(1,-2.9, 1.6, 0.5,  -1,
                    1,-0.9, -1,   -3,   2,
                    1, 1.9, 3.9,  .3,  -2),3,5,byrow=TRUE)
  x1<-runif(100,0,2.5)
  z1<-rnorm(100,0,2.5)
  dat <- data.frame(y=rep(1,100),x1,x2=x1^2,z1,z2=z1^2)
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
                                              dat,theta,dist="negative_binomial")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])

  ## first lets test the classic bernoulli model.
  ## test that silly inputs return the right error message.
  # first test the formulas
  # sam_form_wrong <- y ~ 1 + x1 + x2 + z1 + z2
  # sp_form_wrong <- ~ 1 + x1
  # testthat::expect_error(fm1 <- species_mix(sam_form_wrong, sp_form, 'a', distribution = 'bernoulli', n_mixtures=3))
  # testthat::expect_error(fm2 <- species_mix(sam_form, sp_form_wrong, model_data, distribution = 'bernoulli', n_mixtures=3))
  #
  # # test that null formula fails.
  # testthat::expect_error(fmnb <- species_mix(sam_form, sp_form, model_data, distribution = 'negative_binomial', n_mixtures=3, control = species_mix.control(quiet=TRUE)))

  # set.seed(42)
  # sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  # theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  # dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  # dat[,-1] <- scale(dat[,-1])
  # simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="negative_binomial")
  # model_data <- make_mixture_data(species_data = simulated_data$species_data,
  #                                 covariate_data = simulated_data$covariate_data[,-1])
  #
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

  # test_nb <- nbglm::glm.fit.nbinom(X,y[,3],offset,weights,est_var = FALSE)
  # test_nb <- nbglm::glm.fit.nbinom(X,y[,1],offset,weights,est_var = TRUE)


  # # test a single negative_binomial model
  i <- 1
  testthat::expect_length(ecomix:::apply_glm_sam_inits(i, y, X, site_spp_weights, offset, y_is_na, disty),3)
  fm_negative_binomialint <- surveillance::plapply(1:S, ecomix:::apply_glm_sam_inits,  y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  testthat::expect_length(do.call(cbind,fm_negative_binomialint)[1,],S)

  fm_negative_binomialint
  # test that the starting values work.
  testthat::expect_length(tmp <- ecomix:::fitmix_EM_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control),8)


  # #get the taus
  # starting_values <- ecomix:::initiate_fit_sam(y, X, site_spp_weights, offset, y_is_na, G, S, disty, control)
  # fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,disp=starting_values$disp)
  # first_fit <- list(x = X, y = y, weights=site_spp_weights, offset=offset)
  #
  # # get the loglikelihood based on these values
  # logls <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  # pis <- rep(1/G, G)
  # taus <- ecomix:::get_taus(pis, logls, G, S)
  # taus <- ecomix:::skrink_taus(taus, max_tau=1/G + 0.1, G)
  #
  # ## get to this in a bit
  # gg <- 1
  # testthat::expect_length(ecomix:::apply_glm_group_tau_sam(gg, y, X, site_spp_weights, offset, y_is_na, disty, taus),3)
  #
  # # ## now let's try and fit the optimisation
  # start_vals <- ecomix:::get_starting_values_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
  #
  # tmp <- ecomix:::sam_optimise(y,X,offset,spp_weights,site_spp_weights, y_is_na, S, G, nrow(y), disty, start_vals, control)
  # testthat::expect_length(tmp,15)
  #
  #
  # set.seed(42)
  # sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  # sp_form <- ~1
  # theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  # dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  # dat[,-1] <- scale(dat[,-1])
  # simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="negative_binomial")
  # model_data <- make_mixture_data(species_data = simulated_data$species_data,
  #                                 covariate_data = simulated_data$covariate_data[,-1])
  # fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'negative_binomial',
  #                    n_mixtures=4)
  # testthat::expect_s3_class(fm1,'negative_binomial')
  # testthat::expect_s3_class(fm1,'species_mix')
  #
  # fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'negative_binomial',
  #                    n_mixtures=4,control=species_mix.control(em_prefit = FALSE))
  # testthat::expect_s3_class(fm2,'negative_binomial')
  # testthat::expect_s3_class(fm2,'species_mix')
})
