context('species_mix generic functions four: gaussian functions')
library(ecomix)

testthat::test_that('species mix gaussian', {

  library(ecomix)
  rm(list=ls())
  set.seed(42)
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
  alpha <- runif(50,10,100)
  beta <- matrix(c(3.6,0.5,-0.9,1,4.9,2.9,0.2,-0.4),4,2,byrow=TRUE)
  dat <- data.frame(x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
  simulated_data <- species_mix.simulate(archetype_formula = sam_form,
                                         species_formula = ~1, data = dat,
                                         nArchetypes = 4, alpha = alpha,
                                         beta=beta, family  = "gaussian")
  y <- as.matrix(simulated_data[,grep("spp",colnames(simulated_data))])
  colSums(y>0)
  X <- simulated_data[,-grep("spp",colnames(simulated_data))]
  U <- NULL
  # U <- X[,4,drop=FALSE]
  # X <- X[,-4, drop=FALSE]
  W <- as.data.frame(X[,1,drop=FALSE])
  X <- as.data.frame(X[,-1])

  offset <- rep(0,nrow(y))
  weights <- rep(1,nrow(y))

  spp_weights <- rep(1,ncol(y))
  site_spp_weights <- matrix(1,nrow(y),ncol(y))
  y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  G <- 4
  S <- ncol(y)
  control <- species_mix.control()
  disty <- 6
  size <- rep(1,nrow(y))
  powers <- attr(simulated_data,"powers") # yeah baby
  options(warn=1)
  fm <- ecomix:::get_initial_values_sam(y, X, W, U, site_spp_weights, offset, y_is_na, G, S, disty, size, powers, control = species_mix.control())

  ## now let's try and fit the optimisation
  start_vals <- ecomix:::starting_values_wrapper(y, as.data.frame(X), as.data.frame(W), U, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, size, powers, control)

  tmp <- ecomix:::sam_optimise(y, X, W, U, offset, spp_weights, site_spp_weights, y_is_na, S, G, disty, size, powers, start_vals = start_vals, control)
  testthat::expect_length(tmp,20)

  set.seed(123)
  tmp <- ecomix:::species_mix.fit(y=y, X=as.data.frame(X), W=as.data.frame(W), U=U, G=G, S=S,
                                  spp_weights=spp_weights,
                                  site_spp_weights=site_spp_weights,
                                  offset=offset, disty=disty, y_is_na=y_is_na, size=size, powers=powers,
                                  control=species_mix.control(print_cpp_start_vals = TRUE))

  sp_form <- ~1
  fm1 <- species_mix(archetype_formula = sam_form, species_formula = sp_form,
                     data = simulated_data, family = 'gaussian',
                     nArchetypes = 4)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, data = simulated_data, family = 'gaussian',
                     nArchetypes = 4,control=species_mix.control(ecm_prefit = FALSE),
                     standardise = FALSE)
  testthat::expect_s3_class(fm2,'species_mix')
})
