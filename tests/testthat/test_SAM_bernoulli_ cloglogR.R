context('species_mix generic functions two: bernoulli cloglog functions')
library(ecomix)

testthat::test_that('species mix cloglog bernoulli', {

  library(ecomix)
  rm(list = ls())
  set.seed(42)
  sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:50),
                                                      collapse = ','),
                                       ")~1+x1+x2"))
  sp_form <- ~ 1
  beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,7.9),3,2,byrow=TRUE)
  dat <- data.frame(y=rep(1,200),x1=stats::runif(200,0,2.5),
                    x2=stats::rnorm(200,0,2.5))
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- species_mix.simulate(archetype_formula=sam_form,
                                         species_formula=sp_form,
                                         data = dat,
                                         beta=beta,
                                         family=binomial(link = 'logit'))
  sp_form <- ~1
  fm1 <- species_mix(archetype_formula = sam_form, species_formula = sp_form,
                     data = simulated_data, family = binomial(link = "logit"),
                     nArchetypes = 3)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, data = simulated_data, family = binomial(link = "cloglog"),
                     nArchetypes = 3,control=list(em_prefit = FALSE))
  testthat::expect_s3_class(fm2,'species_mix')

})

