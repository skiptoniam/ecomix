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
  fm1 <- ecomix:::apply_glm_sam_inits(ss, y, X, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,3)

  fm_bern <- surveillance::plapply(seq_len(S), ecomix:::apply_glm_sam_inits, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)

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
  fm1 <- ecomix:::apply_glm_sam_inits(ss, y, X, site_spp_weights, offset, y_is_na, disty)
  testthat::expect_is(fm1,'list')
  testthat::expect_length(fm1,3)

  fm_pois <- surveillance::plapply(seq_len(S), ecomix:::apply_glm_sam_inits, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet)

  alphas <- lapply(fm_pois, `[[`, 1)
  testthat::expect_length(unlist(alphas),S)

  betas <- lapply(fm_pois, `[[`, 2)
  testthat::expect_length(do.call(rbind, betas),S*nP)

  disp <- unlist(lapply(fm_pois, `[[`, 3))
  testthat::expect_true(all(is.na(disp)))


})

testthat::test_that('testing species_mix call', {

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
  testthat::expect_message(fm1 <- species_mix(NULL, sp_form, model_data, distribution = 'bernoulli',
                     n_mixtures=3))

  dup_spp_data <- cbind('spp1'=model_data[,1],model_data)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',c(1,1:20)),collapse = ','),")~1+x1+x2"))

  testthat::expect_message(fm1 <- species_mix(sam_form, sp_form, dup_spp_data,
                                              distribution = 'bernoulli',
                                              n_mixtures=3))
})

testthat::test_that('testing species_mix.multifit call', {

  library(ecomix)
  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))
  sp_form <- ~ 1
  theta <- matrix(c(1,-2.9,-3.6,1,-0.9,1,1,.9,1.9),3,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
                                              dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fmods <- species_mix.multifit(archetype_formula = sam_form,
                                species_formula = sp_form,
                                data = model_data,
                                distribution = 'bernoulli',
                                nstart = 10, n_mixtures=3)
  testthat::expect_is(fmods,'list')
  testthat::expect_length(fmods,10)
  testthat::expect_s3_class(fmods[[1]],'species_mix')

  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))

  sp_form <- ~ 1
  theta <- matrix(c(1,-2.9,-3.6,1,-0.9,1,1,.9,1.9),3,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
                                              dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  testthat::expect_message(fm1 <- species_mix.multifit(NULL, sp_form, model_data, distribution = 'bernoulli',
                                              n_mixtures=3))

  dup_spp_data <- cbind('spp1'=model_data[,1],model_data)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',c(1,1:20)),collapse = ','),")~1+x1+x2"))

  testthat::expect_message(fm1 <- species_mix.multifit(sam_form, sp_form, dup_spp_data,
                                              distribution = 'bernoulli',
                                              n_mixtures=3))

})

testthat::test_that('testing species mix S3 class functions', {

  library(ecomix)
  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))

  sp_form <- ~ 1
  theta <- matrix(c(1,-2.9,-3.6,1,-0.9,1,1,.9,1.9),3,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
                                              dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
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
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm1 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'bernoulli', n_mixtures=4)

  vcv_mat <- vcov(object = fm1)
  testthat::expect_equal(nrow(vcv_mat),nrow(vcv_mat))
  testthat::expect_is(vcv_mat,'matrix')
  testthat::expect_true(all(is.finite(sqrt(diag(vcv_mat)))))

  vcv_mat_bb <- vcov(object = fm1,method = 'BayesBoot', nboot = 100)
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
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="bernoulli")
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
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="poisson")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm2 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'poisson', n_mixtures=4)

  preds3 <- predict(fm2)
  testthat::expect_length(preds3,2)
  testthat::expect_is(preds3,'list')

  preds4 <- predict(fm2, newobs = dat2)
  testthat::expect_is(preds4,'list')

  # negative binomial
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="negative_binomial")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm3 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'negative_binomial', n_mixtures=4)

  preds5 <- predict(fm3)
  testthat::expect_length(preds3,2)
  testthat::expect_is(preds5,'list')

  preds6 <- predict(fm2, newobs = dat2)
  testthat::expect_is(preds6,'list')

  #gaussian
  simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="gaussian")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])
  fm4 <- species_mix(sam_form, species_formula = ~1, model_data, distribution = 'gaussian', n_mixtures=4)

  preds7 <- predict(fm4)
  testthat::expect_length(preds7,2)
  testthat::expect_is(preds7,'list')

  preds8 <- predict(fm4, newobs = dat2)
  testthat::expect_is(preds8,'list')

  testthat::expect_error(preds8 <- predict(fm4, newobs = data.frame(1,dat2)))

})


# testthat::test_that('species mix estimate groups', {
#
#
#   library(ecomix)
#   set.seed(42)
#   sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:100),collapse = ','),")~1+x1+x2"))
#   theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
#   dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
#   dat[,-1] <- scale(dat[,-1])
#   simulated_data <- simulate_species_mix_data(sam_form,~1,dat,theta,dist="bernoulli")
#   model_data <- make_mixture_data(species_data = simulated_data$species_data,
#                                   covariate_data = simulated_data$covariate_data[,-1])
#
#   fm_grps <- species_mix_estimate_groups(archetype_formula = sam_form, species_formula = ~1, data = model_data, distribution = 'bernoulli', n_mixtures=2:10)
#
#
#
#
#
#
#
# })

