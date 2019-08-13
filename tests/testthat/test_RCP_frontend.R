context('regional_mix front end testing, what goes into RCP?')
library(ecomix)

testthat::test_that('testing regional mix front end testing, what goes into RCP?', {

  rm(list = ls())
  set.seed( 151)
  n <- 100
  S <- 10
  nRCP <- 3
  my.dist <- "bernoulli"
  X <- as.data.frame( cbind( x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10)))
  Offy <- log( runif( n, min=30, max=60))
  pols <- list()
  pols[[1]] <- poly( X$x1, degree=3)
  #important to scale covariates so that regimix can get half-way decent starting values
  pols[[2]] <- poly( X$x2, degree=3)
  X <- as.matrix( cbind( 1, X, pols[[1]], pols[[2]]))
  colnames( X) <- c("const", 'x1', 'x2', paste( "x1",1:3,sep='.'), paste( "x2",1:3,sep='.'))
  p.x <- ncol( X[,-(2:3)])
  p.w <- 3
  W <- matrix(sample( c(0,1), size=(n*p.w), replace=TRUE), nrow=n, ncol=p.w)
  colnames( W) <- paste( "w",1:3,sep=".")
  alpha <- rnorm( S,0)
  tau.var <- 0.5
  b <- sqrt( tau.var/2)
  # a double exponential for RCP effects
  tau <- matrix( rexp( n=(nRCP-1)*S, rate=1/b) - rexp( n=(nRCP-1)*S, rate=1/b), nrow=nRCP-1, ncol=S)
  beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2, -76.0, -14.2, -28.3,
                          -36.8, -17.8, -92.9,-2.7), nrow=nRCP-1, ncol=p.x)
  gamma <- matrix( rnorm( S*p.w), ncol=p.w, nrow=S)

  simDat <- regional_mix.simulate(nRCP=nRCP, S=S, p.x=p.x,p.w = p.w, n=n,
                                       alpha=alpha, tau=tau, beta=beta, gamma=gamma,
                                       X=X[,-(2:3)], W=W, distribution=my.dist,
                                       offset=Offy)

  my.form.RCP <- paste( paste( paste(
    'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
    ')',sep=''),
    '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
  my.form.spp <- ~w.1+w.2+w.3
  testthat::expect_null(regional_mix(rcp_formula = NULL, species_formula = my.form.spp, data = simDat,
                                   distribution = "bernoulli", nRCP = 3, inits = unlist( fm.clean[[goodUn]]$coef),
                                   control=list(optimise=FALSE), offset=offset))

  #test two spp have the same name.
  my.form.RCP <- paste( paste( paste(
  'cbind(', paste( paste( 'spp', c(1,1:S), sep=''), collapse=','), sep=''),
  ')',sep=''),
  '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')

   testthat::expect_null(regional_mix(rcp_formula = my.form.RCP,
                                      species_formula = my.form.spp, data = simDat,
                                   distribution = "bernoulli", nRCP = 3,
                                   offset=offset))

  # test one rcp
  my.form.RCP <- paste( paste( paste(
  'cbind(', paste( paste( 'spp', c(1:S), sep=''), collapse=','), sep=''),
  ')',sep=''),
  '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
  regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp, data = simDat,
             distribution = "bernoulli", nRCP = 1, offset=offset)


})


