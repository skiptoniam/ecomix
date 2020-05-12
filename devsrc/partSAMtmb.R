## partial SAM gam
library(ecomix)
set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
spp_form <- as.formula(~1+w1+w2)
beta <- matrix(c(-3.6,0.5,
                 -0.9,1.0,
                 0.9,-2.9,
                 2.2,5.4),
               4,2,byrow=TRUE)
gamma <- matrix(c(rnorm(50,1),rnorm(50,-2)),50,2)
dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- species_mix.simulate(sam_form, spp_form, dat = dat,
                                       beta = beta, gamma = gamma,
                                       n_mixtures = 4,
                                       distribution = "negative_binomial")

archetype_formula <- sam_form
species_formula <- spp_form

test_dat <- ecomix:::clean_data_sam(simulated_data, archetype_formula, species_formula, distribution = 'negative_binomial')
test_dat$mf.X
test_dat$mf.W

y <- simulated_data[,1:50]
X <- ecomix:::get_X_sam(archetype_formula, test_dat$mf.X)
W <- ecomix:::get_W_sam(species_formula, test_dat$mf.W)

n_mixtures <- G <- 4
S <- ncol(y)
spp_weights <- rep(1,S)
site_spp_weights <- matrix(1,nrow(y),S)
disty <- 4
y_is_na <- is.na(y)
inits <- NULL
control <- species_mix.control(quiet = FALSE)
offset <- rep(0,nrow(X))

library(mgcv)
gam_form <- y~s(x0)+s(x1)+s(x2)+s(x3)
archetype_formula <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~s(x1,bs='cs')+s(x2,bs='cs')"))
species_formula <- as.formula(~1+s(w1,bs='cs')+s(w2,bs='cs'))

forms <- list('archetype_formula'=archetype_formula,'species_formula'=species_formula)
dat <- simulated_data

mm <- makeMatrices(forms,dat)

makeMatrices <- function(forms, dat) {

  # results list
  res <- list()
  a <- model.frame(update.formula(forms[["archetype_formula"]], . ~ 1), data = dat)
  Y <- model.extract(a, "response")
  res$Y <- Y
  # remove intercept from archetype formula and whack in a dummy y
  sam.form <- update.formula(forms[["archetype_formula"]], y ~ -1 + .)
  spp.form <- update.formula(forms[["species_formula"]],y ~ 1 +. )

  ## whack in a dummy response variable
  dat$y <- dat[,1]
  gam_sam <- gam(sam.form, data = dat, method = "REML", fit = TRUE)
  res$gam_sam <- gam_sam
  # accumulate smoothing matrix, block diagonal
  if(length(gam_sam$smooth) > 0){
    S_sam_list <- list()
    for(i in seq_along(gam_sam$smooth)){
      smoo <- gam_sam$smooth[[i]]
      for(j in seq_along(smoo$S)){
        S_sam_list <- c(S_sam_list, as(smoo$S[[j]], "sparseMatrix"))
        res$S_sam_n <- c(res$S_sam_n, nrow(smoo$S[[j]]))
        res$s_sam_k <- c(res$s_sam_k, smoo$bs.dim)
        res$s_sam_names <- c(res$s_sam_names, attr(smoo$sp, "names"))
      }
    }
    # build a block diagonal matrix
    res$S_sam <- Matrix::bdiag(S_sam_list)
  }

  samdat <- dat
  ## build design matrix
  mm <- predict(gam_sam, newdata = samdat, type = "lpmatrix")
  # fixed effects design matrix
  if(gam_sam$nsdf>0){
    # fixed effects design matrix
    res$Xs_sam <- mm[, 1:gam_sam$nsdf, drop=FALSE]
    # smooth design matrix
    res$A_sam <- mm[, -c(1:gam_sam$nsdf), drop=FALSE]
  } else {
    res$Xs_sam <- NULL
    res$A_sam <- mm
  }

  # species matricies
  # GAM setup
  gam_spp <- gam(spp.form, data = dat, method = "REML")
  res$gam_spp <- gam_spp
  # accumulate smoothing matrix, block diagonal
  if(length(gam_spp$smooth) > 0){
    S_spp_list <- list()
    for(i in seq_along(gam_spp$smooth)){
      smoo <- gam_spp$smooth[[i]]
      for(j in seq_along(smoo$S)){
        S_spp_list <- c(S_spp_list, as(smoo$S[[j]], "sparseMatrix"))
        res$S_spp_n <- c(res$S_spp_n, nrow(smoo$S[[j]]))
        res$s_spp_k <- c(res$s_spp_k, smoo$bs.dim)
        res$s_spp_names <- c(res$s_spp_names, attr(smoo$sp, "names"))
      }
    }
    # build a block diagonal matrix
    res$S_spp <- bdiag(S_spp_list)
  }

  sppdat <- dat

  ## build design matrix
  mm <- predict(gam_surface, newdata = sppdat, type = "lpmatrix")
  if(gam_spp$nsdf>0){
    # fixed effects design matrix
    res$Xs_spp <- mm[, 1:gam_spp$nsdf, drop=FALSE]
    # smooth design matrix
    res$A_spp <- mm[, -c(1:gam_spp$nsdf), drop=FALSE]
  } else {
    res$Xs_spp <- NULL
    res$A_spp <- mm
  }

  return(res)
}

data <- simulated_data

partialGAMSAM <- function(archetype_formula, species_formula, data, n_mixtures = 3,
                          inits = "random", vcov=FALSE, theta.range=c(0.001,10),
                          pen.parm=1.25, iters=c(10,200), update.kappa=c(1,0.5,1),
                          contr=list(eps=1e-8,init.sd=1), init.fit = NULL){

  oldWarn <- options()$warn
  options(warn = -1)

  require( splines)
  require( MASS)
  #	require( glm2)
  require( mgcv)

  #get data object
  # allData <- getDataObjs1( form.noMix, form.Mix, data)
  forms <- list('archetype_formula'=archetype_formula,'species_formula'=species_formula)
  mm <- makeMatrices(forms,data)

  ###setting up some data structures etc
  call <- match.call()
  G <- as.integer(n_mixtures)
  S <- ncol(mm$Y)
  N <- nrow(mm$Y)

  cat("\nData (and model) summary. There are:\n\t",N," observations,\n\t",S," species, which are grouped into,\n\t",G," species-archetypes. There are\n\t",ncol(mm$A_spp)," smooth parameters for each species, \n\t",ncol(mm$Xs_spp)," fixed effects for species, and \n\t", ncol(mm$A_sam)," smooth parameters for each archetype\n\n")


  ## upto here:
  #Get initial model fits and groupings
  if( is.null( init.fit))
    fits <- initiate.mod.nb( allData, form.noMix, form.Mix, inits, G, theta.range, contr$init.sd)
  else{
    fits <- init.fit
    cat( "Using initial fit object without checking -- good luck\n")
  }

  #the first E-step
  #get the initial pis == taus #well kinda
  pi <- rep( 1/G, G)	#summary( as.factor( fits$grps)) / length( fits$grps)
  logls <- get.logls.nb( allData, fits, G)
  taus <- get.taus( pi, logls, G, S)
  taus <- shrink.taus( taus, maxTau=1/G + 0.1, G)

  logl.old <- logl.new <- -.Machine$double.xmax
  kount <- 1

  cat( "iteration: 0\n")
  print( fits$coefs)
  maxit <- iters[2]
  n.init.steps <- iters[1]

  while( (kount <= maxit & !converged(logl.old, logl.new, eps=contr$eps)) | kount <= n.init.steps){
    old.fits <- fits
    pi <- colSums( taus) / S
    #update the mixing params
    tmp <- nlminb( start=fits$coefs, objective=incom.logl2.nb, gradient=NULL, hessian=NULL, allData=allData, pis=pi, fits=fits, G=G, S=S)
    fits$coefs <- update.coefs( fits$coefs, tmp$par, kapp=update.kappa[3])
    #update the spp-specific params
    out1 <- kronecker( rep( 1, G), allData$outcomes)
    X1.noMix <- kronecker( rep( 1, G), allData$X.noMix)
    offy1 <- kronecker( rep( 1, G), allData$offset)
    offy2 <- allData$X.Mix %*% t( fits$coefs)
    offy2 <- as.numeric( offy2)
    for( ss in 1:S){
      Mix.taus <- rep( taus[ss,], each=N)
      tmpform <- as.formula( paste( paste( 'out1', '[,ss]', sep=''), '-1+X1.noMix+offset( offy1)+offset( offy2)', sep='~'))
      #			fm <- try( glm2( tmpform, weights=Mix.taus, family=negative.binomial( theta=fits$sppTheta[ss], link='log')))
      fm <- try( gam( tmpform, weights=Mix.taus, family=negbin(theta=fits$sppTheta[ss])))
      kount1 <- 1
      while( class( fm) %in% 'try-error' & kount1 < 10){
        kount1 <- kount1 + 1
        theta <- 10 * fits$sppTheta[ss]
        #	  		fm <- try( glm2( tmpform, weights=Mix.taus, family=poisson( link='log')))
        fm <- try( gam( tmpform, weights=Mix.taus, family=negbin(theta=fits$sppTheta[ss])))
      }
      fits$sppCoefs[ss,] <- update.coefs( fits$sppCoefs[ss,], fm$coef, kappa=update.kappa[1])
    }
    #update the dispersion params
    if( kount >= n.init.steps){
      for( ss in 1:S){
        thet <- optimise(f = theta.logl, interval = theta.range, allData=allData, fits=fits, pi=pi, ss=ss, theta.range=theta.range, maximum=TRUE)$maximum
        fits$sppTheta[ss] <- update.coefs( fits$sppTheta[ss], thet, kappa=update.kappa[2])
      }
    }
    cat( "iteration: ", kount, "\n")
    print( fits$coefs)
    #E-step
    logls <- get.logls.nb( allData, fits, G)	###############need to sort out dispersions!##########
    taus <- get.taus( pi, logls, G, S)

    logl.old <- logl.new
    logl.new <- get.incomplete.logl.nb( pi, allData, fits, G, S, theta.range, pen.parm)
    cat( logl.new, '\n')
    kount2 <- 1
    max.decrease <- 100
    cat( "\n")
    kount <- kount + 1
  }


  fits$tau <- taus
  if( kount > maxit)
    cat( "NOT ")
  cat( "converged\n")

  npTot <- (G-1) +  #pi
    prod( dim( fits$sppCoefs)) + #spp parameters
    S +  #spp dispersion
    prod( dim( fits$coefs))  #grp coefs
  BIC <- -2*logl.new + npTot * log( S)

  flag <- TRUE
  if( vcov){
    if( kount < maxit){
      flag <- FALSE
      cat( "Finding covariance matrix for estimates\n")
      allPars <- c( addLogit(pi), as.double( fits$sppCoefs), fits$sppTheta, as.double( fits$coefs))
      tmp1 <- incom.logl2( allPars, allData, G, S, theta.range, pen.parm=pen.parm)
      #      dyn.load( "partialSAMnb.so")
      tmp2 <- .Call( "calcDerHess", as.numeric( allData$outcomes), as.numeric( allData$offset), as.numeric( allPars), as.numeric( allData$X.noMix), as.numeric( allData$X.Mix), as.integer(S), as.integer(G), as.integer(nrow( allData$outcomes)), as.integer( ncol( allData$X.noMix)), as.integer( ncol( allData$X.Mix)), as.numeric( theta.range), as.numeric( pen.parm))
      tmpLogl <- tmp2[1]
      tmpScores <- tmp2[1+1:length( allPars)]
      tmpHess <- matrix( tmp2[1+length( allPars) + 1:(length(allPars)^2)], nrow=length( allPars), ncol=length( allPars))
      try( vcov <- solve( -tmpHess), silent=TRUE)
    }
  }
  if( flag){
    cat( "Not finding covariance matrix for estimates.  Either unrequested or maximum iterations reached.\n")
    vcov <- tmpHess <- NULL
  }

  #	print( fits)

  print( "Done!")
  options( warn=as.numeric( oldWarn))
  return( list( MixCoefs=fits$coefs, sppCoefs=fits$sppCoefs, sppTheta=fits$sppTheta, allData=allData, logl=logl.new, taus=taus, pi=pi, maxit=maxit, niters=kount-1, vcov=vcov, Hessian=tmpHess, BIC=BIC, fits=fits))
  #  return( list( fits=fits, allData=allData))
}
