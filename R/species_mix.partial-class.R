### Main species mix partial functions to export ###

species_mix.partial <- function(archetype_formula, species_formula, data,
                                n_mixtures = 3, offset = NULL, weights = NULL,
                                bb_weights = NULL, control = NULL,
                                inits = "random", standardise = FALSE,
                                titbits = TRUE, #vcov=FALSE,
                                theta.range=c(0.001,10), pen.parm=1.25,
                                iters=c(10,200), update.kappa=c(1,0.5,1),
                                contr=list(eps=1e-8,init.sd=1), init.fit = NULL){

  data <- as.data.frame(data)
  control <- set_control_sam(control)
  if(!control$quiet)
    message( "SAM modelling")
  call <- match.call()
  if(!is.null(archetype_formula)){
    archetype_formula <- stats::as.formula(archetype_formula)
  } else{
    if(!control$quiet)
      message("There is no SAM model! Please provide a model (intercept at least) -- exitting now")
    return(NULL)
  }
  if(!is.null(species_formula))
    species_formula <- stats::as.formula(species_formula)

  mf <- match.call(expand.dots = FALSE)
  if(distribution=="ippm"){
    m <- match(c("data","offset"), names(mf), 0L)
  } else {
    m <- match(c("data","offset","weights"), names(mf), 0L)
  }

  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))
  dat <- clean_data_sam(mf, archetype_formula, NULL, distribution)

  #get data object

  dat <- ecomix:::clean_data_sam(mf, archetype_formula, species_formula, distribution)
  allData <- get_partial_data_objects(species_formula, archetype_formula, dat)

  ###setting up some data structures et
  call <- match.call()
  G <- as.integer(n_mixtures)
  S <- allData$S
  N <- allData$N

  cat("\nData (and model) summary. There are:\n\t",N," observations,\n\t",S," species, which are grouped into,\n\t",G," species-archetypes. There are\n\t",ncol(allData$X.noMix)," parameters for each species, and\n\t",ncol(allData$X.Mix)," parameters for each archetype\n\n")

  #Get initial model fits and groupings
  if( is.null( init.fit))
    fits <- initiate.mod.nb( allData, species_formula, archetype_formula, inits, G, theta.range, contr$init.sd)
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

  # flag <- TRUE
  # if( vcov){
  #   if( kount < maxit){
  #     flag <- FALSE
  #     cat( "Finding covariance matrix for estimates\n")
  #     allPars <- c( addLogit(pi), as.double( fits$sppCoefs), fits$sppTheta, as.double( fits$coefs))
  #     tmp1 <- incom.logl2( allPars, allData, G, S, theta.range, pen.parm=pen.parm)
  #     #      dyn.load( "partialSAMnb.so")
  #     tmp2 <- .Call( "calcDerHess", as.numeric( allData$outcomes), as.numeric( allData$offset), as.numeric( allPars), as.numeric( allData$X.noMix), as.numeric( allData$X.Mix), as.integer(S), as.integer(G), as.integer(nrow( allData$outcomes)), as.integer( ncol( allData$X.noMix)), as.integer( ncol( allData$X.Mix)), as.numeric( theta.range), as.numeric( pen.parm))
  #     tmpLogl <- tmp2[1]
  #     tmpScores <- tmp2[1+1:length( allPars)]
  #     tmpHess <- matrix( tmp2[1+length( allPars) + 1:(length(allPars)^2)], nrow=length( allPars), ncol=length( allPars))
  #     try( vcov <- solve( -tmpHess), silent=TRUE)
  #   }
  # }
  # if( flag){
  #   cat( "Not finding covariance matrix for estimates.  Either unrequested or maximum iterations reached.\n")
  #   vcov <- tmpHess <- NULL
  # }

  #	print( fits)

  print( "Done!")
  options( warn=as.numeric( oldWarn))
  return( list( MixCoefs=fits$coefs, sppCoefs=fits$sppCoefs, sppTheta=fits$sppTheta, allData=allData, logl=logl.new, taus=taus, pi=pi, maxit=maxit, niters=kount-1, vcov=vcov, Hessian=tmpHess, BIC=BIC, fits=fits))
  #  return( list( fits=fits, allData=allData))
}


get_partial_data_objects <- function( f.noMix, f.Mix, dat){

  #organising the data into the two pieces (for spp-specific GLMs and for group-specific GLMS)
  #this is a bit tedious.
  #keeping the names from the two bits
  mix.var.names <- colnames(get_all_vars( f.Mix, dat))
  nomix.var.names <- colnames( get_all_vars( f.noMix, dat))
  #hacking together a combined formula
  tmp.form <- as.formula( paste( f.noMix[2], paste( f.noMix[3], f.Mix[2], sep=" + "), sep=' ~ '))
  #getting combined model.frame
  mf <- model.frame( form=tmp.form, data=dat, na.action=na.pass)
  notNA.pattern <- !apply( mf, 1, function(x) any( is.na( x)))
  mf <- mf[notNA.pattern,]
  #getting the outcomes
  outcomes <- model.response(mf)
  #getting the offset
  offy <- model.offset( mf)
  #getting design matrices -- tedious, I don't know why I can't do this from the existing mf -- it just doesn't work.
  tmp.form1 <- tmp.form; tmp.form1[[2]] <- tmp.form1[[3]]; tmp.form1[[3]] <- NULL
  tmp <- model.matrix( tmp.form1, dat[notNA.pattern,])
  mix.id <- NA
  for( ii in mix.var.names)
    mix.id <- c( mix.id, grep( ii, colnames( tmp)))
  mix.id <- mix.id[!is.na( mix.id)]
  X.Mix <- tmp[,mix.id,drop=FALSE]
  nomix.id <- NA
  for( ii in nomix.var.names)
    nomix.id <- c( nomix.id, grep( ii, colnames( tmp)))
  nomix.id <- nomix.id[!is.na( nomix.id)]
  X.noMix <- tmp[,nomix.id,drop=FALSE]
  if( "(Intercept)" %in% colnames( tmp)){
    tmp.names <- colnames( X.noMix)
    X.noMix <- cbind( 1, X.noMix)
    colnames( X.noMix) <- c( "(Intercept)", tmp.names)
  }
  #complete data frame
  allData <- get_all_vars( tmp.form, dat)
  allData <- allData[notNA.pattern,]

  res <- list( outcomes=outcomes, offset=offy, X.noMix=X.noMix, X.Mix=X.Mix, allData=allData, mix.var.names=mix.var.names, nomix.var.names=nomix.var.names, S=ncol( outcomes), N=nrow( outcomes))

  return( res)

}
