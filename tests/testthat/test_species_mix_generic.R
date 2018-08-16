context('species_mix generic functions')
library(ecomix)


testthat::test_that('species mix generic', {


  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="bernoulli")

  y <- simulated_data$species_data
  X <- simulated_data$covariate_data
  offset <- rep(0,nrow(y))
  # weights <- rep(1,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- length(simulated_data$pi)
  S <- length(simulated_data$sp.int)
  nP <- ncol(X[,-1])
  control <- species_mix.control()

  # test a new glmnet function bernoulli
  ss <- 1
  disty <- 1
  fm1 <- ecomix:::apply_glmnet_sam(ss, y, X, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,3)

  fm_bern <- surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)

  alphas <- lapply(fm_bern, `[[`, 1)
  testthat::expect_length(unlist(alphas),S)

  betas <- lapply(fm_bern, `[[`, 2)
  testthat::expect_length(do.call(rbind, betas),S*nP)

  disp <- unlist(lapply(fm_bern, `[[`, 3))
  testthat::expect_true(all(is.na(disp)))

  ## poisson
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="poisson")

  y <- simulated_data$species_data
  X <- simulated_data$covariate_data
  offset <- rep(0,nrow(y))
  # weights <- rep(1,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- length(simulated_data$pi)
  S <- length(simulated_data$sp.int)
  nP <- ncol(X[,-1])
  control <- species_mix.control()

  ss <- 1
  disty <- 2
  fm1 <- ecomix:::apply_glmnet_sam(ss, y, X, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,3)

  fm_pois <- surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)

  alphas <- lapply(fm_pois, `[[`, 1)
  testthat::expect_length(unlist(alphas),S)

  betas <- lapply(fm_pois, `[[`, 2)
  testthat::expect_length(do.call(rbind, betas),S*nP)

  disp <- unlist(lapply(fm_pois, `[[`, 3))
  testthat::expect_true(all(is.na(disp)))

  ##ippm test
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="ippm")
  y <- simulated_data$species_data
  X <- simulated_data$covariate_data
  offset <- simulated_data$offset
  site_spp_weights <- simulated_data$background_weights
  y_is_na <- simulated_data$y_is_na

  # dim(y)
  # dim(X)
  # dim(y_is_na)
  # dim(site_spp_weights)
  # length(offset)

  ss <- 1
  disty <- 3
  fm1 <- ecomix:::apply_glmnet_sam(ss, y, X, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,3)

  fm_ippm <- surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)

  alphas <- lapply(fm_ippm, `[[`, 1)
  testthat::expect_length(unlist(alphas),S)

  betas <- lapply(fm_ippm, `[[`, 2)
  testthat::expect_length(do.call(rbind, betas),S*nP)

  disp <- unlist(lapply(fm_ippm, `[[`, 3))
  testthat::expect_true(all(is.na(disp)))


  ## negative binomial
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="negative_binomial")

  y <- simulated_data$species_data
  X <- simulated_data$covariate_data
  offset <- rep(0,nrow(y))
  # weights <- rep(1,nrow(y))
  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- length(simulated_data$pi)
  S <- length(simulated_data$sp.int)
  nP <- ncol(X[,-1])
  control <- species_mix.control()

  ss <- 1
  disty <- 4
  fm1 <- ecomix:::apply_glmnet_sam(ss, y, X, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,3)
  #
  fm_nb <- surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)
  #
  alpha <- lapply(fm_nb, `[[`, 1)
  testthat::expect_length(unlist(alpha),S)

  beta <- lapply(fm_nb, `[[`, 2)
  testthat::expect_length(do.call(rbind, beta),S*nP)

  disp <- unlist(lapply(fm_nb, `[[`, 3))
  testthat::expect_length(disp,S)

})


testthat::test_that('species mix generic vcov functions', {

  # build and test a single model.
  # estimate variance-covariance matrix
  library(ecomix)
  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))

  sp_form <- ~ 1
  theta <- matrix(c(1,-2.9,-3.6,1,-0.9,1,1,.9,7.9),3,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
                                              dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli',
   n_mixtures=3)

  vcv_mat <- vcov(object = fm1)
  testthat::expect_equal(nrow(vcv_mat),nrow(vcv_mat))
  testthat::expect_is(vcv_mat,'matrix')
  testthat::expect_true(all(is.finite(sqrt(diag(vcv_mat)))))

  vcv_mat_bb <- vcov(object = fm1,method = 'BayesBoot', nboot = 100)
  testthat::expect_equal(nrow(vcv_mat),nrow(vcv_mat))
  testthat::expect_is(vcv_mat_bb,'matrix')
  testthat::expect_true(all(is.finite(sqrt(diag(vcv_mat_bb)))))

})
