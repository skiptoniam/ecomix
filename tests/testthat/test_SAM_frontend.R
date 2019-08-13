context('species_mix generic functions one')
library(ecomix)

testthat::test_that('species mix internal functions classes work', {


  set.seed(42)
  archetype_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~x1+x2"))
  beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,7.9),3,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(archetype_formula = archetype_form,
                                         species_formula = ~1,
                                         dat = dat,
                                         beta=beta,
                                         n_mixtures = 3,
                                         dist="bernoulli")
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

})
