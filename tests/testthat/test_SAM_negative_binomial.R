context('species_mix generic functions negative binomial functions')



testthat::test_that('species_mix negative binomial', {

  library(ecomix)
  set.seed(42)
  nsp <- 100
  sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:nsp),collapse = ','),")~x1+x2"))
  alpha <- rnorm(nsp,-0.5, .5)
  beta <- matrix(c(3.6,-3.6,
                   -2.5,-2.5,
                   1,-4.5),
                 3,2,byrow=TRUE)
  # delta <- -.4
  x <- runif(200,-2.5,2.5)
  # u <- rnorm(400)
  xpred <- seq(-2.5,2.5,length.out = 100)
  upred <- matrix(seq(-2.5,2.5,length.out = 100),ncol=1)
  matplot(xpred,(exp((cbind(xpred,xpred^2)%*%t(beta)))))
  matplot(xpred,(exp((cbind(xpred,xpred^2)%*%t(beta))-c(upred%*%delta))))
  # matlines(xpred,)

  dat <- data.frame(y=1, x1=x, x2=I(x)^2)#, u)
  simulated_data <- species_mix.simulate(archetype_formula = sam_form,
                                         species_formula = ~1,
                                         all_formula = NULL,
                                         dat = dat,
                                         nArchetypes = 3,
                                         alpha=alpha,
                                         beta=beta,
                                         # delta = delta,
                                         family = "negative.binomial")
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
  G <- 3
  S <- ncol(y)
  control <- species_mix.control()
  disty <- 4
  size <- rep(1,nrow(y))
  powers <- attr(simulated_data,"powers") # yeah baby
  options(warn=1)
  fm <- ecomix:::fit.ecm.sam(y, X, W, U, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, size, powers, control = species_mix.control(em_refit = 3,em_steps = 100))

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
                     data = simulated_data, family = 'negative.binomial',
                     nArchetypes = 3)
  testthat::expect_s3_class(fm1,'species_mix')

  fm2 <- species_mix(sam_form, sp_form, data = simulated_data, family = 'negative.binomial',
                     nArchetypes = 3,control=species_mix.control(em_prefit = FALSE),
                     standardise = FALSE)
  testthat::expect_s3_class(fm2,'species_mix')
})

