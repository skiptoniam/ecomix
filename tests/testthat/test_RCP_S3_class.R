context('regional_mix bernoulli functions')
library(ecomix)

testthat::test_that('testing regional mix S3 functions', {

  library(ecomix)
  rm(list = ls())
  set.seed( 151)
  n <- 100
  S <- 10
  nRCP <- 3
  my.dist <- "negative.binomial"
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
  alpha <- rnorm( S)
  tau.var <- 0.5
  b <- sqrt( tau.var/2)
  #a double exponential for RCP effects
  tau <- matrix( rexp( n=(nRCP-1)*S, rate=1/b) - rexp( n=(nRCP-1)*S, rate=1/b), nrow=nRCP-1, ncol=S)
  beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2, -76.0, -14.2, -28.3,
                          -36.8, -17.8, -92.9,-2.7), nrow=nRCP-1, ncol=p.x)
  gamma <- matrix( rnorm( S*p.w), ncol=p.w, nrow=S)
  logDisp <- log( rexp( S, 1))
  set.seed(121)
  simDat <- regional_mix.simulate( nRCP=nRCP, S=S, p.x=p.x, p.w=p.w, n=n,
                                   alpha=alpha, tau=tau, beta=beta, gamma=gamma,
                                   X=X[,-(2:3)], W=W, family=my.dist,
                                   logDisp=logDisp, offset=Offy)

  #fit the model
  my.form.RCP <- paste( paste( paste(
    'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
    ')',sep=''),
    '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
  my.form.spp <- ~w.1+w.2+w.3
  fm2 <- regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp,
                     data = simDat, family =  "negative.binomial", nRCP = 3, inits = "random2")

  testthat::expect_is(AIC(fm2),'numeric')

  testthat::expect_is(BIC(fm2),'numeric')

  cdres <- cooks.distance(fm2,times=10)
  testthat::expect_s3_class(cdres,"regiCooksD")

  plot(fm2)
  plot(fm2,fitted.scale = "log")
  plot(fm2, type = "deviance" )
  plot(fm2, type = "deviance" , alpha.conf = c(0.75))
  plot(fm2, type = "deviance" , species = fm2$names$spp[1])

  print(fm2)

  resres <- residuals(fm2)
  testthat::expect_is(resres,'matrix')

})


testthat::test_that('Test the prediction functions in RCP',{

  library(ecomix)
  rm(list = ls())
  set.seed( 151)
  n <- 100
  S <- 10
  nRCP <- 3
  my.dist <- "negative.binomial"
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
  alpha <- rnorm( S)
  tau.var <- 0.5
  b <- sqrt( tau.var/2)
  #a double exponential for RCP effects
  tau <- matrix( rexp( n=(nRCP-1)*S, rate=1/b) - rexp( n=(nRCP-1)*S, rate=1/b), nrow=nRCP-1, ncol=S)
  beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2, -76.0, -14.2, -28.3,
                          -36.8, -17.8, -92.9,-2.7), nrow=nRCP-1, ncol=p.x)
  gamma <- matrix( rnorm( S*p.w), ncol=p.w, nrow=S)
  logDisp <- log( rexp( S, 1))
  set.seed(121)
  simDat <- regional_mix.simulate( nRCP=nRCP, S=S, p.x=p.x, p.w=p.w, n=n,
                                   alpha=alpha, tau=tau, beta=beta, gamma=gamma,
                                   X=X[,-(2:3)], W=W, family=my.dist,
                                   logDisp=logDisp, offset=Offy)

  #fit the model
  my.form.RCP <- paste( paste( paste(
    'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
    ')',sep=''),
    '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
  my.form.spp <- ~w.1+w.2+w.3
  fm2 <- regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp,
                     data = simDat, family =  "negative.binomial", nRCP = 3, inits = "random2")
  preds <- predict(fm2)
  testthat::expect_is(preds,'matrix')
  bootmod <- regional_mix_boot(fm2,nboot = 10)
  bootmod1 <- regional_mix_boot(fm2,nboot = 10,type = "SimpleBoot")
  predsboot <- predict(fm2, bootmod)
  predsboot1 <- predict(fm2, bootmod1)
  testthat::expect_is(predsboot,'list')
  testthat::expect_is(predsboot1,'list')

  xwnew <- rbind(replicate(2,cbind(X,W))[,,1],replicate(2,cbind(X,W))[,,2])
  preds_new <- predict(fm2,newdata = xwnew)

  extractAIC(fm2)

  fm2 <- regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp,
                      data = simDat, family =  "negative.binomial",
                      nRCP = 3, inits = "random2",
                      control = list(getScores=TRUE))
})


testthat::test_that('stability',{

  library(ecomix)
  rm(list = ls())
  set.seed( 151)
  n <- 100
  S <- 10
  nRCP <- 3
  my.dist <- "negative.binomial"
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
  alpha <- rnorm( S)
  tau.var <- 0.5
  b <- sqrt( tau.var/2)
  #a double exponential for RCP effects
  tau <- matrix( rexp( n=(nRCP-1)*S, rate=1/b) - rexp( n=(nRCP-1)*S, rate=1/b), nrow=nRCP-1, ncol=S)
  beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2, -76.0, -14.2, -28.3,
                          -36.8, -17.8, -92.9,-2.7), nrow=nRCP-1, ncol=p.x)
  gamma <- matrix( rnorm( S*p.w), ncol=p.w, nrow=S)
  logDisp <- log( rexp( S, 1))
  set.seed(121)
  simDat <- regional_mix.simulate( nRCP=nRCP, S=S, p.x=p.x, p.w=p.w, n=n,
                                   alpha=alpha, tau=tau, beta=beta, gamma=gamma,
                                   X=X[,-(2:3)], W=W, family=my.dist,
                                   logDisp=logDisp, offset=Offy)

  #fit the model
  my.form.RCP <- paste( paste( paste(
    'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
    ')',sep=''),
    '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
  my.form.spp <- ~w.1+w.2+w.3
  fm2 <- regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp,
                     data = simDat, family =  "negative.binomial",
                     nRCP = 3, inits = "random2")
  stab <- stability.regional_mix(fm2,oosSizeRange=c(1,5,10), mc.cores=1, times = 10, doPlot=FALSE)
  plot(stab)

  testthat::expect_error(residuals(fm2,type='RQRS.sim'))
  # test_res <- residuals(fm2,type='RQR.sim',mc.cores=1)
  # testthat::expect_is(test_res,"matrix")

  vcov(fm2)
  tmp <- vcov(fm2,method="FiniteDifference")
  vcov(fm2,method="BayesBoot",nboot=10)
  vcov(fm2,method="SimpleBoot",nboot=10)
  vcov(fm2,method="EmpiricalInfo")

  fm2$vcov <- tmp
  # regional_mix_bootParametric(fm2, nboot = 10)
  summary(fm2)

})


testthat::test_that("RCP membership",{

  library(ecomix)
  rm(list = ls())
  set.seed( 151)
  n <- 100
  S <- 10
  nRCP <- 3
  my.dist <- "poisson"
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
  alpha <- rnorm( S)
  tau.var <- 0.5
  b <- sqrt( tau.var/2)
  #a double exponential for RCP effects
  tau <- matrix( rexp( n=(nRCP-1)*S, rate=1/b) - rexp( n=(nRCP-1)*S, rate=1/b), nrow=nRCP-1, ncol=S)
  beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2, -76.0, -14.2, -28.3,
                          -36.8, -17.8, -92.9,-2.7), nrow=nRCP-1, ncol=p.x)
  gamma <- matrix( rnorm( S*p.w), ncol=p.w, nrow=S)
  # logDisp <- log( rexp( S, 1))
  set.seed(121)
  simDat1 <- regional_mix.simulate( nRCP=nRCP, S=S, p.x=p.x, n=n, #p.w=p.w,
                                   alpha=alpha, tau=tau, beta=beta, #gamma=gamma,
                                   X=X[,-(2:3)], #W=W,
                                   family=my.dist, offset=Offy)

  simDat2 <- regional_mix.simulate( nRCP=nRCP, S=S, p.x=p.x, p.w=p.w, n=n,
                                    alpha=alpha, tau=tau, beta=beta, gamma=gamma,
                                    X=X[,-(2:3)], W=W,
                                    family=my.dist, offset=Offy)

  #fit the model
  my.form.RCP <- paste( paste( paste(
    'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
    ')',sep=''),
    '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
  my.form.spp <- ~1#w.1+w.2+w.3
  fm1 <- regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp,
                      data = simDat1, family =  "poisson",
                      nRCP = 3, inits = "random2")
  my.form.spp <- ~w.1+w.2+w.3


  my.form.RCP <- paste( paste( paste(
    'cbind(', paste( paste( 'spp', c(1,1:S), sep=''), collapse=','), sep=''),
    ')',sep=''),
    '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
  testthat::expect_null(fm2 <- regional_mix.multifit(rcp_formula = my.form.RCP, species_formula = my.form.spp,
                      data = simDat2, family =  "poisson",
                      nRCP = 3, inits = "random"))
  testthat::expect_null(fm2 <- regional_mix.multifit(rcp_formula = NULL, species_formula = my.form.spp,
                                                     data = simDat2, family =  "poisson",
                                                     nRCP = 3, inits = "random"))

  fm2 <- regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp,
                      data = simDat2, family =  "poisson",
                      nRCP = 3, inits = "random2")

  tmpboot1 <- regional_mix_boot(fm1,nboot = 10)
  regional_mix.species_profile(fm1)
  regional_mix.species_profile(fm1,tmpboot1)

  AIC(fm1,k=NULL)
  extractAIC(fm1,k=NULL)

  predict(fm1)
  plot(fm1)
  plot(fm1,species = fm1$names$spp[1])
  testthat::expect_error(plot(fm1,type='blah'))
  testthat::expect_error(plot(fm1,species='blah'))

})
