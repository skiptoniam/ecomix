### Main species mix partial functions to export ###
#' @name species_mix_partial
#' @rdname species_mix_partial
#' @export

"species_mix_partial" <- function(archetype_formula, species_formula, data,
                                  n_mixtures = 3, distribution="negative_binomial",
                                  offset = NULL, weights = NULL,
                                  bb_weights = NULL, control = NULL,
                                  inits=NULL, standardise = FALSE,
                                  titbits = TRUE){

  data <- as.data.frame(data)
  control <- ecomix:::set_control_sam(control)
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

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  if(distribution=="ippm"){
    mf$na.action <- "na.pass"
  } else {
    mf$na.action <- "na.exclude"
  }

  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())


  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))

  #get data object
  dat <- ecomix:::clean_data_sam(mf, archetype_formula, species_formula, distribution)

  # get responses
  y <- stats::model.response(dat$mf.X)

  # logical matirx needed for removing NAs from response and weights mainly used for ippms.
  y_is_na <- is.na(y)

  # check names of reponses
  S <- ecomix:::check_reponse_sam(y)

  # what is the X matrix (archetype covariates)
  X <- get_X_part_sam(archetype_formula = archetype_formula, mf.X = dat$mf.X)

  # what is the W matrix (species covariates)
  W <- get_W_part_sam(species_formula = species_formula, mf.W = dat$mf.W)

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

  tmp$dist <- disty_cases[disty]

  if(n_mixtures==1){
    tmp$pis <- tmp$pis
  } else {
    tmp$pis <- additive_logistic(tmp$eta)
  }

  #calc posterior porbs and pis.
  if(n_mixtures>1)
    tmp$taus <- calc_post_probs_sam(tmp$pis,tmp$loglloglikeSG)

  tmp$pis <- colSums(tmp$taus)/S

  #Information criteria
  tmp <- calc_info_crit_part_sam(tmp)

  #titbits object, if wanted/needed.
  tmp$titbits <- get_titbits_partial_sam(titbits, y, X, W, spp_weights, site_spp_weights,
                                 offset, y_is_na, archetype_formula,
                                 species_formula, control, disty_cases[disty],
                                 tmp$removed_species)
  class(tmp) <- c("species_mix_partial")
}

"species_mix_partial.fit" <- function(y, X, W, G, S, spp_weights,
                                      site_spp_weights, offset, y_is_na, disty,
                                      control, inits=NULL){

  ## currently only use the ECM. We could move the newton-rapts in the future.
  tmp <- fitmix_ECM_partial_sam(y, X, W, spp_weights, site_spp_weights,
                          offset, y_is_na, G, S, disty, control)
  if(G==1)tmp <- clean_ECM_output_one_group(tmp, G, S, disty)
  return(tmp)

}



"fitmix_ECM_partial_sam" <- function(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){

  ite <- 1
  restart_ite <- 1
  logl_old <- -99999999
  logl_new <- -88888888

  # get starting values
  starting_values <- get_initial_values_partial_sam(y = y, X = X, W = W,
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
  logls_mus <- get_logls_partial_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)

  while(control$em_reltol(logl_new,logl_old) & ite <= control$em_steps){
    if(restart_ite>10){
      message('cannot find good starting values with initialisation and random starting values - please check the number of groups and coefs.')
      break
    }

    pis <- colMeans(taus)

    if (any(pis < sqrt(.Machine$double.eps))) {
      if(restart_ite==1){
        cat('Pis have gone to zero - restarting with new initialisation \n')
        starting_values <- get_initial_values_partial_sam(y = y, X = X, W=W,
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
        fits$alpha <- rnorm(S)
        fits$beta <- matrix(rnorm(G*(ncol(X)-1)),G,(ncol(X)-1))
        fits$gamma <- matrix(rnorm(S*(ncol(W)),S,ncol(W)))
        if (disty%in%c(4,6)) fits$theta <- rep(0.05,S)
      }
      ite <- 1
      restart_ite <- restart_ite + 1
    }

    # m-step
    fm_species_coefs <- surveillance::plapply(seq_len(S),
                                       apply_glm_spp_coefs_partial_sams,
                                       y, X, W, G, taus, site_spp_weights, offset,
                                       y_is_na, disty, fits,
                                       .parallel = control$cores,
                                       .verbose = FALSE)
    #check weights in this.
    alpha <- unlist(lapply(fm_species_coefs, `[[`, 1))
    fits$alpha <- alpha#update_alpha_coefs(fits$alpha,alpha)
    gamma <- do.call(rbind,lapply(fm_species_coefs, `[[`, 2))
    fits$gamma <- gamma#update_gamma_coefs(fits$gamma,gamma)


    ## update the betas
    fm_mix_coefs <- surveillance::plapply(seq_len(G),
                                        apply_glm_mix_coefs_partial_sams,
                                        y, X, W, site_spp_weights,
                                        offset, y_is_na, disty, taus,
                                        fits, logls_mus$fitted,
                                        .parallel = control$cores,
                                        .verbose = FALSE)

    # update the coefs.
    fmix_coefs_mat <- t(do.call(cbind,fm_mix_coefs))
    fits$beta <- ecomix:::update_mix_coefs(fits$beta, fmix_coefs_mat)


    ## need a function here that updates the dispersion parameter.
    if(disty%in%c(4,6)){
      fm_theta <- surveillance::plapply(seq_len(S),
                                        apply_optimise_partial_sam_theta,
                                        first_fit, fits,
                                        G, disty, pis,
                                       .parallel = control$cores,
                                       .verbose = FALSE)
      thetas <- unlist(lapply(fm_theta, `[[`, 1))
      theta <- log(1/thetas)
      fits$theta <- ecomix:::update_sp_dispersion(fits$theta,theta,0.5)
    }

    # e-step
    # get the log-likes and taus
    logls_mus <- get_logls_partial_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)
    taus <- ecomix:::get_taus(pis, logls_mus$logl_sp, G, S)
    # taus <- shrink_taus(taus,)

    #update the likelihood
    logl_old <- logl_new
    logl_new <- get_incomplete_logl_partial_sam(eta = ecomix:::additive_logistic(pis,inv = TRUE)[-G],
                                        first_fit, fits, spp_weights, G, S, disty)
    if(!control$quiet)cat("Iteration ",ite,"\n")
    if(!control$quiet)cat("Loglike: ", logl_new,"\n")
    if(!control$quiet)cat("Pis: ", pis,"\n")
    ite <- ite + 1
  }

  taus <- data.frame(taus)
  names(taus) <- paste("grp.", seq_len(G), sep = "")
  alpha_out <- fits$alpha
  beta_out <- fits$beta
  gamma_out <- fits$gamma
  theta_out <- fits$theta
  names(pis) <- paste("G", seq_len(G), sep = ".")
  eta <- additive_logistic(pis, TRUE)[-1]

  # estimate log-likelihood
  # logl_new <- get_incomplete_logl_partial_sam(eta, first_fit, fits, spp_weights, G, S, disty)

  return(list(logl = logl_new, alpha = alpha_out, beta = beta_out,
              gamma = gamma_out, theta = theta_out, eta = eta, pis = pis,
              taus = round(taus,4), first_fit = first_fit))

}


"apply_glm_spp_coefs_partial_sams" <- function(ss, y, X, W, G, taus,
                                               site_spp_weights,
                                               offset, y_is_na, disty, fits){

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
    outcomes <- as.numeric(y[ids_i,ss])
  }
  out1 <- kronecker(rep( 1, G), outcomes)
  X1 <- kronecker(rep( 1, G), X[ids_i,])
  W1 <- kronecker(rep( 1, G), cbind(1,W[ids_i,]))
  wts1 <- kronecker(rep( 1, G),
                    as.numeric(site_spp_weights[ids_i,ss]))*rep(taus[ss,],
                                                                each=length(site_spp_weights[ids_i,ss]))
  offy1 <- kronecker(rep( 1, G), offset[ids_i])
  offy2 <- X[ids_i,-1] %*% t(fits$beta)
  offy2 <- as.numeric(offy2)
  offy <- offy1 + offy2

  if(disty %in% c(1,2,3,6)){
    ft_sp <- try(stats::glm.fit(x=as.data.frame(W1),
                                y=as.numeric(out1),
                                weights=as.numeric(wts1),
                                offset=as.numeric(offy),
                                family=fam), silent=FALSE)
    if (class(ft_sp) %in% 'try-error'){
      print(paste0(ss,"\n"))
      my_coefs <- rep(NA, ncol(X1))
    } else {
      my_coefs <- coef(ft_sp)
      names(my_coefs) <- c(colnames(y)[ss],colnames(W))
    }
  }
  if(disty %in% 4){
    ft_sp <- glm.fit.nbinom(x=as.matrix(W1),
                            y=as.numeric(out1),
                            weights=as.numeric(wts1),
                            offset=as.numeric(offy))
    my_coefs <- ft_sp$coef
    names(my_coefs) <- c(colnames(y)[ss],colnames(W))
  }
  return(list(alpha = my_coefs[1], gamma = my_coefs[-1]))
}

"apply_glm_mix_coefs_partial_sams" <- function(gg, y, X, W, site_spp_weights, offset,
                                                 y_is_na, disty,
                                                 taus, fits, mus){

  ### setup the data stucture for this model.
  Y_taus <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
  X_no_NA <- list()
  W_no_NA <- list()
  for (jj in 1:ncol(y)){
    X_no_NA[[jj]] <- X[!y_is_na[,jj],-1,drop=FALSE]
    W_no_NA[[jj]] <- cbind(1,W[!y_is_na[,jj],,drop=FALSE])
  }
  X_taus <- do.call(rbind, X_no_NA)
  W_taus <- do.call(rbind, W_no_NA)
  n_ys <- sapply(X_no_NA,nrow)
  wts_taus <- rep(taus[,gg,drop=FALSE],c(n_ys))
  if(disty==4)wts_taus <- rep(taus[,gg],c(n_ys))/(1+rep(exp(-fits$theta),n_ys)*as.vector(mus[gg,,]))

  site_weights <- as.matrix(as.matrix(unlist(as.data.frame(site_spp_weights[!y_is_na]))))
  wts_tausXsite_weights <- wts_taus*site_weights
  offy_mat <- replicate(ncol(y),offset)
  offy1 <- unlist(as.data.frame(offy_mat[!y_is_na]))
  offy2 <- cbind(1,W) %*% t(cbind(fits$alpha,fits$gamma))
  offy2 <- unlist(as.data.frame(offy2[!y_is_na]))
  offy <- as.numeric(offy1 + offy2)

  # which family to use?
  if( disty == 1)
    fam <- binomial()
  if( disty == 2 | disty == 3 |disty == 4)
    fam <- poisson()
  if( disty == 6)
    fam <- gaussian()
  if (disty==3){
    Y_taus <- as.matrix(Y_taus/site_weights)
  } else {
    Y_taus <- as.matrix(Y_taus)
  }

  if(disty %in% c(1,2,3,4,6)){ #don't use for tweedie
    ft_mix <- try(glm.fit(x = as.data.frame(X_taus),
                          y = as.numeric(Y_taus),
                          weights = c(wts_tausXsite_weights),
                          offset = offy,
                          family = fam), silent = TRUE)
    if (class(ft_mix) %in% 'try-error'){
      mix_coefs <- rep(NA, ncol(X_taus[,-1]))
    } else {
      mix_coefs <- ft_mix$coefficients
    }

  }
  return(c(mix_coefs))
}


"apply_optimise_partial_sam_theta" <- function(ss, first_fit, fits,
                                               G, disty, pis,
                                               theta.range = c(0.0001,100)){
  thet <- optimise(f = theta_logl_part_sam, interval = theta.range, ss, first_fit,
                   fits, G, disty, pis, theta.range,
                   maximum = TRUE)$maximum
  return(thet)
}

"theta_logl_part_sam" <- function(theta, ss, first_fit, fits, G,
                          disty, pis, theta.range){

  eta.species <- cbind(1,first_fit$W) %*% t(cbind(fits$alpha[ss],fits$gamma[ss,,drop=FALSE]))
  eta.mixture <- first_fit$x[,-1] %*% t(fits$beta)
  logls <- rep(NA,G)
  for(gg in seq_len(G)) {
    eta.all <- eta.species + eta.mixture[,gg] + first_fit$offset
    if(disty == 4){
      link <- make.link('log')
      logls[gg] <- sum(dnbinom(first_fit$y[,ss], mu = link.fun$linkinv(eta.all), size = 1/exp(theta), log = TRUE));
    }
    if(disty == 6){
      link <- make.link('identity')
      logls[gg] <- sum(dnorm(first_fit$y[,ss], mean = link.fun$linkinv(eta.all), sd = exp(theta), log = TRUE));
    }
  }

  ak <- logls + log( pis)
  am <- max(ak)
  ak <- exp( ak-am)
  sppLogls <- am + log( sum( ak))

  pen.max <- theta.range[2]
  pen.min <- theta.range[1]
  shape1 <- shape2 <- 1.25

  if(disty == 4) sppLogls <- sppLogls + dbeta( (theta-pen.min) / (pen.max-pen.min), shape1, shape2, log=TRUE)

  return(sppLogls)
}

"get_initial_values_partial_sam" <- function(y, X, W, spp_weights, site_spp_weights,
                                     offset, y_is_na, G, S, disty, control){

  # get intial model fits
  starting_values <- initiate_fit_partial_sam(y, X, W, spp_weights, site_spp_weights,
                                      offset, y_is_na, G, S, disty, control)

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

"get_logls_partial_sam" <- function(first_fit, fits, spp_weights, G, S, disty, get_fitted=TRUE){

  if(get_fitted) fitted_values <- array(0,dim=c(G,nrow(first_fit$y),S))
  # if(is.null(spp_weights))spp_weights <- rep(1,S) #for bayesian boostrap.

  logl_sp <- matrix(NA, nrow=S, ncol=G)
  eta.species <- as.matrix(cbind(1,first_fit$W)) %*% t(as.matrix(cbind(fits$alpha,fits$gamma)))
  eta.mixture <- as.matrix(first_fit$x[,-1]) %*% t(fits$beta)

  #bernoulli
  if(disty == 1){
    link <- make.link('logit')
    for( ss in seq_len(S)){
      for( gg in seq_len(G)){
      eta <- eta.species[,ss] + eta.mixture[,gg] + first_fit$offset
      if(get_fitted)fitted_values[gg,,ss] <- link$linkinv(eta)
      logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(eta), log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  #poisson
  if(disty==2){
    link <- make.link('log')
    for(ss in seq_len(S)){
      for(gg in seq_len(G)){
        eta <- eta.species[,ss] + eta.mixture[,gg] + first_fit$offset
        if(get_fitted)fitted_values[gg,,ss] <- link$linkinv(eta)
        logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss],link$linkinv(eta),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  #ippm
  if(disty==3){
    link <- make.link('log')
    for(ss in seq_len(S)){
      sp_idx<-!first_fit$y_is_na[,ss]
      for(gg in seq_len(G)){
        #lp is the same as log_lambda (linear predictor)
        eta <- cbind(first_fit$x[sp_idx,1],first_fit$W[sp_idx,]) %*% c(fits$alpha[ss],fits$gamma[ss,]) + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$beta[gg,] + first_fit$offset[sp_idx]
        if(get_fitted) fitted_values[gg,sp_idx,ss] <- link$linkinv(eta)
        logl_sp[ss,gg] <- (first_fit$y[sp_idx,ss] %*% eta - first_fit$site_spp_weights[sp_idx,ss] %*% link$linkinv(eta))
      }
    }
  }

  #negative binomial
  if(disty==4){
    link <- make.link('log')
    for(ss in seq_len(S)){
      for(gg in seq_len(G)){
        eta <- eta.species[,ss] + eta.mixture[,gg] + first_fit$offset
        if(get_fitted)fitted_values[gg,,ss] <- link$linkinv(eta)
        logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss], mu=link$linkinv(eta), size=1/exp(-fits$theta[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  # gaussian
  if(disty==6){
    link <- make.link('identity')
    for(ss in seq_len(S)){
      for(gg in seq_len(G)){
        eta <- eta.species[,ss] + eta.mixture[,gg] + first_fit$offset
        if(get_fitted)fitted_values[gg,,ss] <- link$linkinv(eta)
        logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=link$linkinv(eta),sd=exp(fits$theta[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  out.list <- list(logl_sp=logl_sp)
  if(get_fitted) out.list$fitted = fitted_values
  return(out.list)
}


"get_incomplete_logl_partial_sam" <- function(eta, first_fit, fits, spp_weights, G, S, disty){

  pis <- ecomix:::additive_logistic(eta)
  logl_sp <- matrix(NA, nrow=S, ncol=G)
  eta.species <- cbind(1,first_fit$W) %*% t(cbind(fits$alpha,fits$gamma))
  eta.mixture <- first_fit$x[,-1] %*% t(fits$beta)

  #bernoulli
  if(disty == 1){
    link <- make.link('logit')
    for( ss in seq_len(S)){
      for( gg in seq_len(G)){
        eta <- eta.species[,ss] + eta.mixture[,gg] + first_fit$offset
        logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(eta), log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  #poisson
  if(disty==2){
    link <- make.link('log')
    for(ss in seq_len(S)){
      for(gg in seq_len(G)){
        eta <- eta.species[,ss] + eta.mixture[,gg] + first_fit$offset
        logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss],link$linkinv(eta),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  #ippm
  if(disty==3){
    link <- make.link('log')
    for(ss in seq_len(S)){
      sp_idx<-!first_fit$y_is_na[,ss]
      for(gg in seq_len(G)){
        #lp is the same as log_lambda (linear predictor)
        eta <- cbind(first_fit$x[sp_idx,1],first_fit$W[sp_idx,]) %*% c(fits$alpha[ss],fits$gamma[ss,]) + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$beta[gg,] + first_fit$offset[sp_idx]
        logl_sp[ss,gg] <- (first_fit$y[sp_idx,ss] %*% eta - first_fit$site_spp_weights[sp_idx,ss] %*% link$linkinv(eta))
      }
    }
  }

  #negative binomial
  if(disty==4){
    link <- make.link('log')
    for(ss in seq_len(S)){
      for(gg in seq_len(G)){
        eta <- eta.species[,ss] + eta.mixture[,gg] + first_fit$offset
        logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss], mu=link$linkinv(eta), size=1/exp(-fits$theta[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  # gaussian
  if(disty==6){
    link <- make.link('identity')
    for(ss in seq_len(S)){
      for(gg in seq_len(G)){
        eta <- eta.species[,ss] + eta.mixture[,gg] + first_fit$offset
        logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=link$linkinv(eta),sd=exp(fits$theta[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
ak <- logl_sp + matrix(rep(log( pis), each=S), nrow=S, ncol=G)
am <- apply( ak, 1, max)
ak <- exp( ak-am)
sppLogls <- am + log( rowSums( ak))
logl <- sum( sppLogls)

return(logl)
}


"get_partial_data_objects" <- function(archetype_formula, species_formula, dat){

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


"get_X_part_sam" <-function(archetype_formula, mf.X){
  form.X <- archetype_formula
  form.X[[2]] <- NULL
  form.X <- stats::as.formula(form.X)
  X <- stats::model.matrix(form.X, mf.X)
  return( X)
}


"get_W_part_sam" <- function(species_formula, mf.W){
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
"apply_glmnet_sam_inits" <- function(ss, y, X, W, site_spp_weights, offset, y_is_na, disty){

  # which family to use?
  if(disty == 1)
    fam <- "binomial" #glmnet
  # fam <- binomial() #glm
  if(disty == 2 | disty == 3 | disty == 4)
    fam <- "poisson"
  # fam <- poisson()
  if(disty == 6)
    fam <- "gaussian"
  # fam <- gaussian()

  ids_i <- !y_is_na[,ss]

  if (disty==3){
    outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  } else {
    outcomes <- as.matrix(y[ids_i,ss])
  }

  if( disty %in% c(1,2,3,4,6)){
    lambda.seq <- sort( unique( c( seq( from=1/0.001, to=1, length=25), seq( from=1/0.1, to=1, length=10))), decreasing=TRUE)

    ft_sp <- try(glmnet::glmnet(y=outcomes, x=as.matrix(cbind(W[ids_i,-1,drop=FALSE],X[ids_i,,drop=FALSE])),
                                family=fam, offset=offset[ids_i],
                                weights=as.matrix(site_spp_weights[ids_i,ss]),
                                alpha=0,
                                lambda=lambda.seq,
                                standardize=FALSE,
                                intercept=TRUE), silent=FALSE)
    locat.s <- 1/1
    my.coefs <- glmnet::coef.glmnet(ft_sp, s=locat.s)
    if( any( is.na( my.coefs))){  #just in case the model is so badly posed that mild penalisation doesn't work...
      my.coefs <- glmnet::coef.glmnet(ft_sp, s=lambda.seq)
      lastID <- apply( my.coefs, 2, function(x) !any( is.na( x)))
      lastID <- tail( (seq_along( lastID))[lastID], 1)
      my.coefs <- my.coefs[,lastID]
    }
    if (any(class(ft_sp) %in% 'try-error')){
      my_coefs <- rep(NA, ncol(X[ids_i,]))
      names(my_coefs) <- colnames(cbind(W[ids_i,,drop=FALSE],X[ids_i,,drop=FALSE]))
    } else {
      my_coefs <- t(as.matrix(my.coefs))
    }

  ##estimate the starting dispersion parameter.
  theta <- NA
  if( disty == 4){
      tmp <- MASS::theta.mm(outcomes, as.numeric(predict(ft_sp, s=locat.s,
                                                         type="response",
                                                         newx=as.matrix(cbind(W[ids_i,-1,drop=FALSE],X[ids_i,,drop=FALSE])),
                                                         newoffset=offset[ids_i])),
                            weights=as.matrix(site_spp_weights[ids_i,ss]),
                            dfr=length(outcomes), eps=1e-4)
      if( tmp>2)
        tmp <- 2
      theta <- log( 1/tmp)
  }
  if( disty == 6){
  preds <- as.numeric( predict(ft_sp, s=locat.s, type="link",
                                   newx=as.matrix(cbind(W[ids_i,-1,drop=FALSE],X[ids_i,,drop=FALSE])), newoffset=offset[ids_i]))
  theta <- log( sqrt( sum((outcomes - preds)^2)/length(outcomes)))  #should be something like the resid standard
    }
  }
  # species intercpets
  alpha <- my_coefs[1]
  # mixture coefs
  beta <- my_coefs[match(colnames(X), colnames(my_coefs))]
  # species coefs apart from intercept
  gamma <-  my_coefs[match(colnames(W), colnames(my_coefs))]

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

  #get taus as starting values
  if(G==1){
    taus <- matrix(1,nrow=ncol(y), ncol = G)
  } else {
    taus <- matrix(0,nrow=ncol(y), ncol= G)
    if(length(sel_omit_spp)>0){
      for(j in 1:length((seq_len(S))[-sel_omit_spp]))
        taus[(seq_len(S))[-sel_omit_spp][j],fmmvnorm$cluster[j]] <- 1
      taus[sel_omit_spp,] <- matrix(runif(length(sel_omit_spp)*G),length(sel_omit_spp), G)
    } else {
      for(j in seq_len(S))taus[j,fmmvnorm$cluster[j]] <- 1
    }
  }
  taus <- taus/rowSums(taus)
  taus <- ecomix:::shrink_taus(taus, max_tau = 1/G + 0.1, G=G)
  pis <- colMeans(taus)


  res <- list()
  res$grps <- tmp_grp
  res$alpha <- alpha
  res$beta <- grp_coefs
  res$beta_all <- beta
  res$gamma <- gamma
  res$theta <- theta
  res$taus <- taus
  res$pis <- pis
  res$species_to_remove <- species_to_remove

  return( res)
}

"calc_info_crit_part_sam" <-  function(tmp) {
  # k <- length(unlist(tmp$coefs))

  k <- length(tmp$eta) +  #pi
    prod(dim(fits$gamma)) + #spp parameters
    length(tmp$alpha) +  #spp dispersion
    prod( dim( fits$beta))  #grp coefs

  tmp$BIC <- -2 * tmp$logl + log(length(tmp$alpha)) * k
  tmp$AIC <- -2 * tmp$logl + 2 * k
  return( tmp)
}


"update_coefs_part_sam" <- function( old, new, kappa=0.5) {
  if( is.na( old))
    tmp <- new
  else
    tmp <- old + kappa*(new-old)
  return( tmp)
}
