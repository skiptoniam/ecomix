context('species_mix generic data sructures')
library(ecomix)

testthat::test_that('species mix one covariate and one group', {

  rm(list = ls())
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
  beta <- matrix(c(-1.6,-2.5,
                   -0.9,1,
                   2.9,2.9,
                   0.2,-0.4),4,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form, ~1, dat = dat,
                                         beta = beta,
                                         nArchetypes = 4,
                                         family="bernoulli")
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
  sp_form <- ~ 1
  all_form <- NULL
  fm1 <- species_mix(sam_form, sp_form, all_form, simulated_data,
                     family = 'bernoulli',
                     nArchetypes = 4)

  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1"))
  sp_form <- ~ 1
  expect_warning(fm1 <- species_mix(sam_form, sp_form, simulated_data,
                                    family = 'bernoulli',
                                    nArchetypes = 4))

  # test with one group
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
  beta <- matrix(c(-1.6,-2.5,
                   -0.9,1,
                   2.9,2.9,
                   0.2,-0.4),4,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form,~1,dat = dat, nArchetypes = 4,
                                         beta = beta,family = "bernoulli")
  sp_form <- ~ 1
  fm1 <- species_mix(sam_form, sp_form, simulated_data,
                     family = 'bernoulli',
                     nArchetypes = 1)

})

testthat::test_that('species mix test multifit', {

  rm(list = ls())
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~x1+x2"))
  beta <- matrix(c(-1.6,-2.5,
                   -0.9,1,
                   2.9,2.9),3,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(sam_form, ~1, dat = dat,
                                         beta = beta,
                                         nArchetypes = 3,
                                         family="bernoulli")
  sp_form <- ~ 1
  fm1 <- species_mix.multifit(sam_form, sp_form, simulated_data,
                              family = 'bernoulli',
                              nArchetypes = 3)

})

