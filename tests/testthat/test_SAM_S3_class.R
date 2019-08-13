context('species_mix generic functions S3 class functions')
library(ecomix)

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

  coef(fm1)
  print(fm1)
  testthat::expect_error(summary(fm1))
  fm1$vcov <- vcov(fm1,method = 'BayesBoot', nboot = 5)
  testthat::expect_is(summary(fm1),'matrix')
  testthat::expect_length(AIC(fm1),1)
  testthat::expect_length(BIC(fm1),1)
  testthat::expect_is(coef(fm1),'list')
  testthat::expect_is(summary(fm1),'matrix')
  testthat::expect_is(predict(fm1),'matrix')
  testthat::expect_is(residuals(fm1),'matrix')
  })



testthat::test_that('species mix generic vcov functions', {

  # build and test a single model.
  # estimate variance-covariance matrix
  library(ecomix)
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:100),collapse = ','),")~1+x1+x2"))
  beta <- matrix(c(1.6,0.5,-0.9,1,2.9,2.9,0.2,-0.4),4,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form, ~1, dat = dat, n_mixtures = 4,
                                         beta = beta, distribution = "bernoulli")
  fm1 <- species_mix(sam_form, species_formula = ~1, simulated_data,
                     distribution = 'bernoulli', n_mixtures=4)

  vcv_mat_bb <- vcov(object = fm1,method = 'BayesBoot', nboot = 10)
  # testthat::expect_equal(nrow(vcv_mat),nrow(vcv_mat))
  testthat::expect_is(vcv_mat_bb,'matrix')
  testthat::expect_true(all(is.finite(sqrt(diag(vcv_mat_bb)))))

})

testthat::test_that('species mix predict functions', {

  # build and test a single model.
  library(ecomix)
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:100),collapse = ','),")~1+x1+x2"))
  beta <- matrix(c(1.6,0.5,-0.9,1,2.9,2.9,0.2,-0.4),4,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form, ~1, dat = dat, n_mixtures = 4,
                                         beta = beta, distribution = "bernoulli")
  fm1 <- species_mix(sam_form, species_formula = ~1, simulated_data, distribution = 'bernoulli', n_mixtures=4)

  preds <- predict(fm1)
  testthat::expect_length(preds,400)
  testthat::expect_is(preds,'matrix')

  dat2 <- data.frame(x1=runif(100,-2.5,2.5),x2=rnorm(100,-2.5,2.5))
  preds2 <- predict(fm1, newobs = dat2)
  testthat::expect_is(preds2,'matrix')

  # poisson
  simulated_data <- species_mix.simulate(sam_form, ~1, dat = dat, n_mixtures = 4,
                                         beta = beta, distribution = "poisson")
  fm2 <- species_mix(sam_form, species_formula = ~1, simulated_data, distribution = 'poisson', n_mixtures=4)

  preds3 <- predict(fm2)
  testthat::expect_length(preds3,400)
  testthat::expect_is(preds3,'matrix')

  preds4 <- predict(fm2, newobs = dat2)
  testthat::expect_is(preds4,'matrix')

  # negative binomial
  simulated_data <- species_mix.simulate(sam_form, ~1, dat = dat, n_mixtures = 4,
                                         beta = beta, distribution = "negative_binomial")
  fm3 <- species_mix(sam_form, species_formula = ~1, simulated_data, distribution = "negative_binomial", n_mixtures=4)

  preds5 <- predict(fm3)
  testthat::expect_length(preds3,400)
  testthat::expect_is(preds5,'matrix')

  preds6 <- predict(fm3, newobs = dat2)
  testthat::expect_is(preds6,'matrix')

  #gaussian
  simulated_data <- species_mix.simulate(sam_form, ~1, dat = dat, n_mixtures = 4,
                                         beta = beta, distribution = "gaussian")
  fm4 <- species_mix(sam_form, species_formula = ~1, simulated_data,
                     distribution = "gaussian", n_mixtures=4)

  preds7 <- predict(fm4)
  testthat::expect_length(preds7,400)
  testthat::expect_is(preds7,'matrix')

  preds8 <- predict(fm4, newobs = dat2)
  testthat::expect_is(preds8,'matrix')

  testthat::expect_error(preds8 <- predict('a'))

})
