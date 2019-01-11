context('regional_mix generic functions')
library(ecomix)


testthat::test_that('testing regional mix function', {

set.seed( 151)
n <- 100
S <- 10
nRCP <- 3
my.dist <- "negative_binomial"
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
simDat <- simulate_regional_mix_data( nRCP=nRCP, S=S, p.x=p.x, p.w=p.w, n=n, alpha=alpha, tau=tau,
                      beta=beta, gamma=gamma, X=X[,-(2:3)], W=W, dist=my.dist, logDisp=logDisp, offset=Offy)

#fit the model
my.form.RCP <- paste( paste( paste(
  'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
  ')',sep=''),
  '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
my.form.spp <- ~w.1+w.2+w.3
fm <- regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp,
                   data = simDat, distribution =  "negative_binomial", nRCP = 3, inits = "random2")

testthat::expect_s3_class(fm,'regional_mix')

## Not run:
#fit the model using multiple starting values
fm <- regional_mix.multifit(rcp_formula = my.form.RCP, species_formula = my.form.spp, data = simDat,
                        dist="negative_binomial", nRCP = 3, inits = "random2", offset=offset, nstart=10, titbits=FALSE,
                        mc.cores=1)

#sometimes the model 'mis-fits' and one or more of the RCP groups has no sites associated
#with it.  These need to be removed (based on the colSums of the posterior probabilities)
postProbSums <- t( sapply( fm, function(x) colSums( x$postProbs)))
#Identify those models with zero posterior prob classes
allGoodUns <- apply( postProbSums, 1, function(x) all(x!=0))
#subset the fits
fm.clean <- fm[allGoodUns]
#choose the model with the lowest BIC
goodUn <- which.min(sapply(fm.clean, ecomix:::BIC.regional_mix))
#Using the 'best' model, use regimix(qv) again to additional model output needed for other
#functions (e.g. plot.regimix(qv), predict.regimix(qv) and regiboot(qv)). Note that the
#model is not estimated again (see control argument of the following regimix(qv) call.
fm.final <- regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp, data = simDat,
                     dist="negative_binomial", nRCP = 3, inits = unlist( fm.clean[[goodUn]]$coef),
                     control=list(optimise=FALSE), offset=offset)
testthat::expect_s3_class(fm.final,'regional_mix')

# test if null formula
testthat::expect_null(regional_mix(rcp_formula = NULL, species_formula = my.form.spp, data = simDat,
                                   dist="negative_binomial", nRCP = 3, inits = unlist( fm.clean[[goodUn]]$coef),
                                   control=list(optimise=FALSE), offset=offset))

#test two spp have the same name.
my.form.RCP <- paste( paste( paste(
  'cbind(', paste( paste( 'spp', c(1,1:S), sep=''), collapse=','), sep=''),
  ')',sep=''),
  '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')

testthat::expect_null(regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp, data = simDat,
                                   dist="negative_binomial", nRCP = 3, inits = unlist( fm.clean[[goodUn]]$coef),
                                   control=list(optimise=FALSE), offset=offset))

# test one rcp
my.form.RCP <- paste( paste( paste(
  'cbind(', paste( paste( 'spp', c(1:S), sep=''), collapse=','), sep=''),
  ')',sep=''),
  '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
regional_mix(rcp_formula = my.form.RCP, species_formula = my.form.spp, data = simDat,
                                   dist="negative_binomial", nRCP = 1, inits = unlist( fm.clean[[goodUn]]$coef),
                                   control=list(optimise=FALSE), offset=offset)

# test one rcp
my.form.RCP <- paste( paste( paste(
  'cbind(', paste( paste( 'spp', c(1:S), sep=''), collapse=','), sep=''),
  ')',sep=''),
  '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
regional_mix(rcp_formula = my.form.RCP, species_formula = NULL, data = simDat,
             dist="negative_binomial", nRCP = 1, inits = unlist( fm.clean[[goodUn]]$coef),
             control=list(optimise=FALSE), offset=offset)


})


testthat::test_that('testing regional mix diagonstic functions', {
library(ecomix)
ecomix:::cooks.distance.regional_mix(fm.final)

})


