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
                                         data = dat,
                                         beta=beta,
                                         nArchetypes = 3,
                                         family = "bernoulli")
  model_data <- simulated_data

  #test formula error
  testthat::expect_error(fm1 <- species_mix(archetype_formula = NA, species_formula = ~1, data = model_data, family = 'bernoulli', nArchetypes = 3))
  testthat::expect_error(fm3 <- species_mix(archetype_formula = NA, species_formula = ~1, data = model_data, family = 'poisson', nArchetypes=3))
  testthat::expect_error(fm4 <- species_mix(archetype_formula = NA, species_formula = ~1, data = model_data, family = 'ippm', nArchetypes=3))
  testthat::expect_error(fm5 <- species_mix(archetype_formula = NA, species_formula = ~1, data = model_data, family = 'negative_binomial', nArchetypes=3))

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

})
