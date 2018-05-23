context('species_mix-class')

testthat::test_that('species mix functions classes work', {

  library(ecomix)
  set.seed(42)
  form <- as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))
  theta <- matrix(c(1,-2.9,-3.6,1,-0.9,1,1,.9,7.9),3,3,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(form,dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                  covariate_data = simulated_data$covariate_data[,-1])

  #test formula error
  testthat::expect_error(fm1 <- species_mix(NA, model_data, distribution = 'bernoulli', n_mixtures=3))
  testthat::expect_error(fm2 <- species_mix(NA, model_data, distribution = 'bernoulli_sp', n_mixtures=3))
  testthat::expect_error(fm3 <- species_mix(NA, model_data, distribution = 'poisson', n_mixtures=3))
  testthat::expect_error(fm4 <- species_mix(NA, model_data, distribution = 'ippm', n_mixtures=3))
  testthat::expect_error(fm5 <- species_mix(NA, model_data, distribution = 'negative_binomial', n_mixtures=3))

  ## test the internal functions
  ## test to see if the species formula checks are working.

  f1 <- y ~ x
  f2 <- y ~ x + z
  f3 <- y ~ 1
  f4 <- NULL

  testthat::expect_true(check_species_formula(f1)==2)
  testthat::expect_true(check_species_formula(f2)==2)
  testthat::expect_true(check_species_formula(f3)==1)
  testthat::expect_true(check_species_formula(f4)==0)

})
