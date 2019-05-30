### Main species mix partial functions to export ###

"species_mix_partial" <- function(archetype_formula, species_formula, data,
                                n_mixtures = 3, offset = NULL, weights = NULL,
                                bb_weights = NULL, control = NULL,
                                inits = "random", standardise = FALSE,
                                titbits = TRUE, #vcov=FALSE,
                                theta.range=c(0.001,10), pen.parm=1.25,
                                iters=c(10,200), update.kappa=c(1,0.5,1),
                                contr=list(eps=1e-8,init.sd=1), init.fit = NULL){

  data <- as.data.frame(data)
  control <- ecomix:::set_control_sam(control)
  if(!control$quiet)
    message( "Partial SAM modelling")
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

  #get data object
  dat <- ecomix:::clean_data_sam(mf, archetype_formula, species_formula, distribution)
  # allData <- get_partial_data_objects(archetype_formula, species_formula, mf)

  # get responses
  y <- stats::model.response(dat$mf.X)

  # logical matirx needed for removing NAs from response and weights mainly used for ippms.
  y_is_na <- is.na(y)

  # check names of reponses
  S <- ecomix:::check_reponse_sam(y)

  # what is the X matrix (archetype covariates)
  X <- ecomix:::get_X_sam(archetype_formula = archetype_formula, mf.X = dat$mf.X)

  # what is the W matrix (species covariates)
  W <- get_W_sam(species_formula = species_formula, mf.W = dat$mf.W)

  x.means <- NULL
  x.sds <- NULL
  w.means <- NULL
  w.sds <- NULL
  if (standardise == TRUE) {
    stand.X <- standardise.X(X[, -1])
    X <- as.matrix(cbind(1, stand.X$X))
    stand.W <- standardise.X(W)
    W <- as.matrix(stand.X$X)
    x.means <- stand.X$dat.means
    x.sds <- stand.X$dat.sds
    w.means <- stand.W$dat.means
    w.sds <- stand.W$dat.sds
  }

  # summarising data to console
  ecomix:::print_input_sam(y, X, W, S, archetype_formula, species_formula, distribution,
                  quiet=control$quiet)

  #get distribution
  disty_cases <- c("bernoulli","poisson","ippm","negative_binomial","tweedie",
                   "gaussian")
  disty <- ecomix:::get_distribution_sam(disty_cases, distribution)

  # get offsets
  offset <- ecomix:::get_offset_sam(dat$mf.X)

  # get the weights
  species_names <- colnames(y)
  site_spp_weights <- ecomix:::get_site_spp_weights_sam(mf, weights, species_names, distribution)
  spp_weights <- ecomix:::check_spp_weights(bb_weights,S)

  if(distribution=='ippm'){
    if(!all(colnames(y)==colnames(site_spp_weights))){
      stop(cat('When modelling a inhomogeneous poisson point process model,
               \n species data colnames must match weights colnames.\n\n
               Species data colnames from "data" are:\n',colnames(y),'.\n\n
               While the colnames of the weights are:\n',
               colnames(site_spp_weights),'\n'))
    }
    if(any(dim(y)!=dim(site_spp_weights))){
      stop('When modelling a inhomogenous poisson point process model,
           weights needs to have the same dimensions at the
           species data - n_sites x n_species')
    }
  }

  # fit partial species mix model
  tmp <- species_mix_partial.fit(y=y, X=X, W=W, G=n_mixtures, S=S,
                                 spp_weights=spp_weights,
                                 site_spp_weights=site_spp_weights,
                                 offset=offset, disty=disty, y_is_na=y_is_na,
                                 control=control, inits=inits)


  # if( is.null( init.fit))
    # fits <- initiate.mod.nb( allData, species_formula, archetype_formula, inits, G, theta.range, contr$init.sd)
  # else{
    # fits <- init.fit
    # cat( "Using initial fit object without checking -- good luck\n")
  # }

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

"species_mix_partial.fit" <- function(y, X, W, G, S, spp_weights, site_spp_weights,
                              offset, y_is_na, disty, control, inits=NULL){

  if(G==1){
    tmp <- fitmix_ECM_sam(y, X, spp_weights, site_spp_weights,
                          offset, y_is_na, G, S, disty, control)
    tmp <- clean_ECM_output_one_group(tmp, G, S, disty)
    return(tmp)
  }

  if(is.null(inits)){
    starting_values  <-  get_starting_values_sam(y = y, X = X, W = W,
                                                 spp_weights = spp_weights,
                                                 site_spp_weights = site_spp_weights,
                                                 offset = offset,
                                                 y_is_na = y_is_na,
                                                 G = G, S = S,
                                                 disty = disty,
                                                 control = control)

  } else {
    if(!control$quiet)message('Be careful! You are using your own initial starting values to optimise the species_mix model.')
    inits <- setup_inits_sam(inits, S=S, G=G, np=ncol(X[,-1]), disty, return_list = TRUE)
    print(inits)
    starting_values <- inits
  }

  tmp <- sam_optimise(y, X, offset, spp_weights, site_spp_weights,
                      y_is_na, S, G, disty, starting_values, control)

  return(tmp)
}



"fitmix_ECM_partial_sam" <- function(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){

  ite <- 1
  restart_ite <- 1
  logl_old <- -99999999
  logl_new <- -88888888

  # get starting values
  starting_values <- get_initial_values_sam(y = y, X = X,
                                            spp_weights = spp_weights,
                                            site_spp_weights = site_spp_weights,
                                            offset = offset, y_is_na = y_is_na,
                                            G = G, S = S,
                                            disty = disty,
                                            control = control)

  # first e-step
  fits <- starting_values$fits
  taus <- starting_values$taus
  pis <- starting_values$pis
  # cat('start pis', pis,'\n')
  first_fit <- starting_values$first_fit
  logls_mus <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)

  while(control$em_reltol(logl_new,logl_old) & ite <= control$em_steps){
    if(restart_ite>10){
      message('cannot find good starting values with initialisation and random starting values - please check the number of groups and coefs.')
      break
    }

    pis <- colMeans(taus)

    if (any(pis < sqrt(.Machine$double.eps))) {
      if(restart_ite==1){
        cat('Pis have gone to zero - restarting with new initialisation \n')
        starting_values <- get_initial_values_sam(y = y, X = X,
                                                  spp_weights = spp_weights,
                                                  site_spp_weights = site_spp_weights,
                                                  offset = offset, y_is_na = y_is_na,
                                                  G = G, S = S,
                                                  disty = disty,
                                                  control = control)
        pis <- starting_values$pis
        fits <- starting_values$fits
        taus <- starting_values$taus
        first_fit <- starting_values$first_fit
      } else {
        cat('Pis have gone to zero - restarting with random inits \n')
        taus <- matrix(runif(S*G),S); taus <- taus/rowSums(taus);
        pis <- colSums(taus)/S;
        fits$beta <- matrix(rnorm(G*(ncol(X)-1)),G,(ncol(X)-1))
        fits$alpha <- rnorm(S)
        if (disty%in%c(4,6)) fits$disp <- rep(0.05,S)
      }
      ite <- 1
      restart_ite <- restart_ite + 1
    }

    if(control$regularisation){
      alpha_estimater <- ecomix:::apply_glmnet_sam_sp_intercepts
      beta_estimater <- ecomix:::apply_glmnet_group_tau_sam
    } else {
      alpha_estimater <- ecomix:::apply_glm_sam_sp_intercepts
      beta_estimater <- ecomix:::apply_glm_group_tau_sam
    }

    # m-step
    fm_sp_int <- surveillance::plapply(seq_len(S),
                                       alpha_estimater,
                                       y, X, G, taus, site_spp_weights, offset,
                                       y_is_na, disty, fits,
                                       .parallel = control$cores,
                                       .verbose = FALSE)
    #check weights in this.
    alpha <- unlist(lapply(fm_sp_int, `[[`, 1))
    fits$alpha <- update_sp_coefs(fits$alpha,alpha)

    ## update the betas
    fmix_coefs <- surveillance::plapply(seq_len(G),
                                        beta_estimater,
                                        y, X, site_spp_weights,
                                        offset, y_is_na, disty, taus,
                                        fits, logls_mus$fitted,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

    # update the coefs.
    fmix_coefs_mat <- t(do.call(cbind,fmix_coefs))
    fits$beta <- update_mix_coefs(fits$beta, fmix_coefs_mat)


    ## need a function here that updates the dispersion parameter.
    if(disty%in%c(4,6)){
      fm_disp <- surveillance::plapply(1:S, apply_optimise_spp_theta,
                                       first_fit, fits,
                                       G, disty, taus,
                                       .parallel = control$cores,
                                       .verbose = FALSE)
      disp <- unlist(lapply(fm_disp, `[[`, 1))
      disp <- log(1/disp)
      fits$disp <- update_sp_dispersion(fits$disp,disp,0.5)
    }

    # e-step
    # get the log-likes and taus
    logls_mus <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
    taus <- get_taus(pis, logls_mus$logl_sp, G, S)
    # taus <- shrink_taus(taus,)

    #update the likelihood
    logl_old <- logl_new
    logl_new <- get_incomplete_logl_sam(eta = additive_logistic(pis,inv = TRUE)[-G],
                                        first_fit, fits, spp_weights, G, S, disty)
    if(!control$quiet)cat("Iteration ",ite,"\n")
    if(!control$quiet)cat("Loglike: ", logl_new,"\n")
    if(!control$quiet)cat("Pis: ", pis,"\n")
    ite <- ite + 1
  }

  taus <- data.frame(taus)
  names(taus) <- paste("grp.", 1:G, sep = "")
  int_out <- fits$alpha
  fm_out <- fits$beta
  names(pis) <- paste("G", 1:G, sep = ".")
  eta <- additive_logistic(pis, TRUE)[-1]

  # estimate log-likelihood
  logl_new <- get_incomplete_logl_sam(eta, first_fit, fits, spp_weights, G, S, disty)

  return(list(logl = logl_new, alpha = int_out, beta = fm_out, disp = fits$disp,
              eta = eta, pis = pis, taus = round(taus,4), first_fit = first_fit))

}


"get_initial_values_partial_sam" <- function(y, X, W, spp_weights, site_spp_weights,
                                     offset, y_is_na, G, S, disty, control){

  # get intial model fits
  starting_values <- initiate_fit_partial_sam(y, X, W, spp_weights, site_spp_weights,
                                      offset, y_is_na, G, S, disty, control)

  # #if any are errors then remove them from the models for ever.
  # updated_y <- update_species_data_structure(y, y_is_na,
  #                                            spp_weights,
  #                                            site_spp_weights,
  #                                            starting_values$species_to_remove)
  # y <- updated_y[[1]]
  # y_is_na <- updated_y[[2]]
  # spp_weights <- updated_y[[3]]
  # site_spp_weights <- updated_y[[4]]
  # S <- ncol(y)

  fits <- list(alpha=starting_values$alpha,
               beta=starting_values$beta,
               gamma=starting_values$gamma,
               theta=starting_values$theta)
  first_fit <- list(x = X, W=W, y = y, spp_weights = spp_weights,
                    site_spp_weights = site_spp_weights,
                    offset = offset, y_is_na = y_is_na,
                    removed_species = starting_values$species_to_remove)

  res <- list()
  res$fits <- fits
  res$first_fit <- first_fit
  res$taus <- starting_values$taus
  res$pis <- colSums(starting_values$taus)/S
  return(res)
}


get_partial_data_objects <- function(archetype_formula, species_formula, dat){

  #organising the data into the two pieces (for spp-specific GLMs and for group-specific GLMS)
  #this is a bit tedious.
  #keeping the names from the two bits
  mix.var.names <- colnames(get_all_vars( archetype_formula, dat))
  nomix.var.names <- colnames( get_all_vars( species_formula, dat))
  #hacking together a combined formula
  tmp_form <- as.formula( paste( archetype_formula[2], paste( species_formula[2], archetype_formula[3], sep=" + "), sep=' ~ '))
  #getting combined model.frame
  mf <- model.frame( form=tmp_form, data=dat, na.action=na.pass)
  notNA.pattern <- !apply( mf, 1, function(x) any( is.na( x)))
  mf <- mf[notNA.pattern,]
  #getting the outcomes
  outcomes <- model.response(mf)
  #getting the offset
  offy <- model.offset(mf)
  #getting design matrices -- tedious, I don't know why I can't do this from the existing mf -- it just doesn't work.
  tmp_form1 <- tmp_form; tmp_form1[[2]] <- tmp_form1[[3]]; tmp_form1[[3]] <- NULL
  tmp <- model.matrix( tmp_form1, dat[notNA.pattern,])
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
  allData <- get_all_vars( tmp_form, dat)
  allData <- allData[notNA.pattern,]

  res <- list( outcomes=outcomes, offset=offy, X.noMix=X.noMix, X.Mix=X.Mix, allData=allData, mix.var.names=mix.var.names, nomix.var.names=nomix.var.names, S=ncol( outcomes), N=nrow( outcomes))

  return( res)
}




"get_W_sam" <- function(species_formula, mf.W){
  form.W <- species_formula
  if( !is.null( species_formula)){
    if( length( form.W)>2)
      form.W[[2]] <- NULL #get rid of outcomes
    W <- model.matrix( form.W, mf.W)
    tmp.fun <- function(x){ all( x==1)}
    intercepts <- apply( W, 2, tmp.fun)
    W <- W[,!intercepts,drop=FALSE]
  }
  else
    W <- -999999
  return( W)
}


## this function should run a glm for each species and collect the species intercept, species covars, mixture covars and any dispersion parameters.
"apply_glm_partial_sam_inits" <- function(ss, y, X, W, site_spp_weights, offset, y_is_na, disty){

  # which family to use?
  if(disty == 1)
    fam <- binomial()
  if(disty == 2 | disty == 3 | disty == 4)
    fam <- poisson()
  if(disty == 6)
    fam <- gaussian()

  ids_i <- !y_is_na[,ss]

  if (disty==3){
    outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  } else {
    outcomes <- as.matrix(y[ids_i,ss])
  }

  if( disty %in% c(1,2,3,6)){
    ft_sp <- try(stats::glm.fit(x=as.data.frame(cbind(X[ids_i,,drop=FALSE],W[ids_i,,drop=FALSE])),
                                y = outcomes,
                                weights=as.numeric(site_spp_weights[ids_i,ss]),
                                offset=offset[ids_i],
                                family=fam), silent=FALSE)
    if (class(ft_sp) %in% 'try-error'){
      my_coefs <- rep(NA, ncol(X[ids_i,]))
    } else {
      my_coefs <- coef(ft_sp)
    }
  }
  disp <- NA
  if(disty == 4){
    ft_sp <- try(glm.fit.nbinom(x=as.matrix(cbind(X[ids_i,,drop=FALSE],W[ids_i,,drop=FALSE])),
                                y=as.numeric(outcomes),
                                weights=as.numeric(site_spp_weights[ids_i,ss]),
                                offset=offset[ids_i],est_var=FALSE), silent = TRUE)

    if (class(ft_sp) %in% 'try-error'){
      my_coefs <- rep(NA, ncol(X[ids_i,]))
    } else {
      my_coefs <- ft_sp$coef
    }
    tmp <- ft_sp$theta
    if(tmp>2) tmp <- 2
    disp <- log(1/tmp)
  }
  if( disty == 6){
    preds <- ecomix:::predict.glm.fit(ft_sp, cbind(X[ids_i,,drop=FALSE],W[ids_i,,drop=FALSE]), offset[ids_i], disty)
    disp <- log(sqrt(sum((outcomes - preds)^2)/length(outcomes)))  #should be something like the resid standard Deviation.
  }

  # need to clean this up so I get a spp vars, mix vars.
  # species intercpets
  alpha <- my_coefs[1]
  # mixture coefs
  beta <- my_coefs[match(colnames(X[,-1]), names(my_coefs))]
  # species coefs apart from intercept
  gamma <-  my_coefs[match(colnames(W), names(my_coefs))]
  # dispersion parameter
  theta <- disp

  return(list(alpha = alpha, beta = beta, gamma = gamma, theta = theta))
}

"initiate_fit_partial_sam" <- function(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){

  # Obtaining initial values via spp-specific GLMs - could play with a GAM option or penalised glm
  fm_sp_mods <-  surveillance::plapply(seq_len(S), apply_glm_partial_sam_inits, y, X, W,
                                       site_spp_weights, offset, y_is_na, disty,
                                       .parallel = control$cores, .verbose = FALSE)

  ## turn lists into vectors or matrices.
  alpha <- unlist(lapply(fm_sp_mods, `[[`, 1))
  beta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))
  gamma <- do.call(rbind,lapply(fm_sp_mods, `[[`, 3))
  theta <- unlist(lapply(fm_sp_mods, `[[`, 4))

  names(alpha) <- colnames(y)
  names(theta) <- colnames(y)

  ## if there are troublesome species remove them
  species_to_remove <- which(apply(beta, 1, function(x) all(is.na(x))))

  if(length(species_to_remove)>0){
    #update fits
    alpha <- alpha[-species_to_remove]
    beta <- beta[-species_to_remove,,drop=FALSE]
    disp <- disp[-species_to_remove]

    # update y, y_is_na and weights
    updated_y <- update_species_data_structure(y, y_is_na, spp_weights, site_spp_weights, species_to_remove)
    y <- updated_y[[1]]
    y_is_na <- updated_y[[2]]
    site_spp_weights <- updated_y[[3]]
  } else {
    species_to_remove <- NA
  }

  ## work out mixing groups only on the more common species?
  prev_min_sites <- control$minimum_sites_occurrence
  sel_omit_spp <- which(colSums(y>0, na.rm = TRUE) <= prev_min_sites)
  if(length(sel_omit_spp)>0) beta <- beta[-sel_omit_spp,,drop=FALSE]

  if(G==1) control$init_method <- 'kmeans'


  if(control$init_method=='kmeans'){
    # if(!control$quiet)message( "Initial groups by K-means clustering\n")
    fmmvnorm <- stats::kmeans(beta, centers=G, nstart=100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G)
  }

  if(control$init_method=='kmed'){
    # message( "Initial groups parameter estimates by K-medoids\n")
    mrwdist <- kmed::distNumeric(beta, beta, method = "mrw")
    fmmvnorm <- kmed::fastkmed(mrwdist, ncluster = G, iterate = 100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- beta[fmmvnorm$medoid,,drop=FALSE]
    grp_coefs <- matrix(grp_coefs,nrow=G)
  }

  if(control$init_method=='random2' | is.null(tmp_grp)){
    fmmvnorm <- stats::kmeans(beta, centers=G, nstart=100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G)

    random_coefs <- sam_random_inits(alpha, grp_coefs, disp, S, G, mult=0.3)
    alpha <- random_coefs[[1]]
    grp_coefs <- random_coefs[[2]]
    if(disty%in%c(4,6)) disp <- random_coefs[[3]]

  }

  if(ncol(X[,-1,drop=FALSE])==1) {
    names(grp_coefs)[2] <- names(X[,-1,drop=FALSE])[2]
  }else{
    colnames(grp_coefs) <- colnames(X[,-1,drop=FALSE])
  }

  res <- list()
  res$grps <- tmp_grp
  res$alpha <- alpha
  res$beta <- grp_coefs
  res$beta_all <- beta
  res$gamma <- gamma
  res$theta <- theta
  return( res)
}
