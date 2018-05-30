context('species_mix bernoulli')

testthat::test_that('species mix functions classes work', {

  library(ecomix)
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2+z1+z2"))
  sp_form <- NULL
  theta <- matrix(c(1,-2.9,3.6,.5,-1,1,-0.9,1,-3,2,1,-.9,7.9,3,-2),3,5,byrow=TRUE)
  x1<-runif(100,0,2.5)
  z1<-rnorm(100,0,2.5)
  dat <- data.frame(y=rep(1,100),x1,x2=x1^2,z1,z2=z1^2)
  dat[,-1] <- scale(dat[,-1])
  simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
                                               dat,theta,dist="bernoulli")
  model_data <- make_mixture_data(species_data = simulated_data$species_data,
                                   covariate_data = simulated_data$covariate_data[,-1])
  testthat::expect_error(fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli', n_mixtures=NULL))
  # fm2 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli', n_mixtures=3)
  # fm3 <- species_mix(sam_form, ~1, model_data, distribution = 'bernoulli', n_mixtures=3)





})
