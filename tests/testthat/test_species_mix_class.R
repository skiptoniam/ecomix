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
  testthat::expect_error(fm1 <- species_mix(NA, model_data, distribution = 'bernoulli_sp', n_mixtures=3))

})
