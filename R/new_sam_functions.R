## generic functions for species mix.

"lambda_penalisation_fun" <- function(x,lambda,kappa=0.1){ #assumes that x spans to pretty-well the unpenalised estiamtes
  min.effective.penalty <- min( which( abs( x-tail( x, 1)) < 0.01 * abs( tail( x, 1))))    #the first that lambda that gives a coef close to the last lambda's corresponding coef
  min.effective.penalty <- lambda[min.effective.penalty]
  target.penalty <- kappa * min.effective.penalty
  res.pos <- which.min( (lambda-target.penalty)^2)
  res <- x[res.pos]
  return( res)
}

"apply_glmnet_sam" <- function(ss, y, X, site_spp_weights, offset, y_is_na, disty){

  options(warn = -1)
  # which family to use?
  if( disty == 1)
    fam <- "binomial"
  if( disty == 2 | disty == 3 | disty == 4)
    fam <- "poisson"
  if( disty == 6)
    fam <- "gaussian"

  ids_i <- !y_is_na[,ss]

  if (disty==3){ outcomes <- as.matrix(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  } else { outcomes <- as.matrix(y[ids_i,ss])
  }

  #lambdas for penalised glm
  lambda.seq <- sort( unique( c( seq( from=1/0.1, to=1, length=10), seq( from=1/0.1, to=1, length=10),seq(from=0.9, to=10^-2, length=10))), decreasing=TRUE)
  if( disty != 5){ #don't use for tweedie
  ft_sp <- glmnet::glmnet(x=as.matrix(X[ids_i,-1]),
                          y=outcomes,
                          weights=c(site_spp_weights[ids_i,ss]),
                          offset=offset[ids_i],
                          family=fam,
                          alpha=0,
                          lambda = lambda.seq,
                          standardize = FALSE,
                          intercept = TRUE)
  my_coefs <- apply(glmnet::coef.glmnet(ft_sp), 1, lambda_penalisation_fun, lambda.seq)
  }
  disp <- NA
  if( disty == 4){
    locat.s <- lambda.seq[max(which(as.matrix(glmnet::coef.glmnet(ft_sp))==my_coefs,arr.ind = TRUE)[,2])]
    preds <-as.numeric(predict(ft_sp, s=locat.s,
                               type="response",
                               newx=X[ids_i,-1],
                               offset=offset[ids_i]))
    tmp <- MASS::theta.mm(outcomes, preds,
                          weights=c(site_spp_weights[ids_i,ss]),
                          dfr=length(y[ids_i,ss]),
                          eps=1e-4)
    if(tmp>2)
      tmp <- 2
      disp <- log( 1/tmp)
  }
  if( disty == 6){
    preds <-as.numeric(predict(ft_sp, s=locat.s,
                               type="response",
                               newx=X[ids_i,-1],
                               offset=offset[ids_i]))
    disp <- log(sqrt(sum((outcomes - preds)^2)/length(outcomes)))  #should be something like the resid standard Deviation.
  }

  return(list(alpha = my_coefs[1], beta = my_coefs[-1], disp = disp))

}

"apply_glm_group_tau_sam" <- function (gg, y, X, site_spp_weights, offset, y_is_na, disty, tau){

    ### setup the data stucture for this model.
    Y_tau <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
    X_no_NA <- list()
    for (jj in 1:ncol(y)){
      X_no_NA[[jj]] <- X[!y_is_na[,jj],]
    }
    X_tau <- do.call(rbind, X_no_NA)
    n_ys <- sapply(X_no_NA,nrow)
    wts_tau <- rep(tau[,gg],c(n_ys))


    ippm_weights <- as.matrix(as.matrix(unlist(as.data.frame(site_spp_weights[!y_is_na]))))
    Z_tau <- as.matrix(Y_tau/ippm_weights)
    wts_tauXippm_weights <- wts_tau*ippm_weights
    offy_mat <- replicate(ncol(y),offset)
    offy <- unlist(as.data.frame(offy_mat[!y_is_na]))

    options(warn = -1)
    # which family to use?
    if( disty == 1)
      fam <- "binomial"
    if( disty == 2 | disty == 3 | disty == 4)
      fam <- "poisson"
    if( disty == 6)
      fam <- "gaussian"

    if (disty==3){ Y_tau <- as.matrix(Y_tau/ippm_weights)
    } else { Y_tau <- as.matrix(Y_tau)
    }

    #lambdas for penalised glm
    lambda.seq <- sort( unique( c( seq( from=1/0.1, to=1, length=10), seq( from=1/0.1, to=1, length=10),seq(from=0.9, to=10^-2, length=10))), decreasing=TRUE)
    if( disty != 5){ #don't use for tweedie
      ft_mix <- glmnet::glmnet(x=as.matrix(X_tau[,-1]),
                               y=as.matrix(Y_tau),
                               weights=c(wts_tauXippm_weights),
                               offset = offy,
                               family=fam,
                               alpha=0,
                               lambda = lambda.seq,
                               standardize = FALSE,
                               intercept = TRUE)
      my_coefs <- apply(glmnet::coef.glmnet(ft_mix), 1, lambda_penalisation_fun, lambda.seq)
    }
    # disp <- NA
    # if( disty == 4){
    #   locat.s <- lambda.seq[max(which(as.matrix(glmnet::coef.glmnet(ft_mix))==my_coefs,arr.ind = TRUE)[,2])]
    #   preds <-as.numeric(predict(ft_mix, s=locat.s,
    #                              type="response",
    #                              newx=X_tau[,-1],
    #                              offset=offy))
    #   tmp <- MASS::theta.mm(Y_tau, preds,
    #                         weights=wts_tauXippm_weights,
    #                         dfr=nrow(Y_tau),
    #                         eps=1e-4)
    #   if(tmp>2)
    #     tmp <- 2
    #   disp <- log( 1/tmp)
    # }
    # if( disty == 6){
    #   preds <-as.numeric(predict(ft_mix, s=locat.s,
    #                              type="response",
    #                              newx=X_tau[,-1],
    #                              offset=offu))
    #   disp <- log(sqrt(sum((Y_tau - preds)^2)/length(Y_tau)))  #should be something like the resid standard Deviation.
    # }

    return(as.matrix(my_coefs))
  }


"get_starting_values_sam" <- function(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){

  temp.warn <- getOption( "warn")
  options( warn=-1)

  #if emfit is in the control do an EM fit to get good starting values for c++
  if(disty==3)control$em_prefit<-FALSE

  if(isTRUE(control$em_prefit)){

    if(!control$quiet)message('Using EM algorithm to find starting values; using',
                              control$em_refit,'refits\n')

    emfits <- list()
    for(ii in seq_len(control$em_refit)){
        emfits[[ii]] <-  fitmix_EM_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
    }
    bf <- which.max(vapply(emfits,function(x)c(x$logl),c(logl=0)))
    emfit <- emfits[[bf]]
    start_vals <- list(alphas=emfit$alpha,
                       betas=emfit$beta,
                       disp=emfit$disp,
                       pis=emfit$pis)
  } else {
    if(!control$quiet)message('You are not using the EM algorith to find starting values; starting values are
                              generated using',control$init_method,'\n')
    starting_values <- get_initial_values_sam(y = y, X = X,
                                              spp_weights = spp_weights,
                                              site_spp_weights = site_spp_weights,
                                              offset = offset, y_is_na = y_is_na,
                                              G = G, S = S,
                                              disty=disty,
                                              control = control)
    start_vals <- list(alphas=starting_values$fits$alphas,
                       betas=starting_values$fits$betas,
                       disp=starting_values$fits$disp,
                       pis=starting_values$pis)
  }

  ## all the things we need to c++ optimisation.
  start_vals$eta <- additive_logistic(start_vals$pis,inv = TRUE)[-G]
  start_vals$nS <- S
  start_vals$nG <- G
  start_vals$nObs <- nrow(y)
  return(start_vals)
}


"get_initial_values_sam" <- function(y, X, spp_weights = NULL, site_spp_weights, offset, y_is_na, G, S, disty, control){

  starting_values <- ecomix:::initiate_fit_sam(y, X, site_spp_weights, offset, y_is_na, G, S, disty, control)
  fits <- list(alphas=starting_values$alphas,betas=starting_values$betas,disp=starting_values$disp)
  first_fit <- list(x = X, y = y, weights = weights, offset = offset, y_is_na = y_is_na)
  logls <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- get_taus(pis, logls, G, S)
  taus <- skrink_taus(taus,max_tau=0.9, G) #max_tau=1/G + 0.1

  #use glm for the step.
  fmix_coefs <- surveillance::plapply(1:G, apply_glm_group_tau_sam,
                                      first_fit$y,
                                      first_fit$x,
                                      first_fit$weights,
                                      first_fit$offset,
                                      first_fit$y_is_na,
                                      disty,
                                      taus,
                                      .parallel = control$cores,
                                      .verbose = !control$quiet)

  #update the mix coefs.
  fmix_coefs <- t(do.call(cbind,fmix_coefs)[-1,])
  fits$betas <- update_mix_coefs(fits$betas,fmix_coefs)

  res <- list()
  res$fits <- fits
  res$first_fit <- first_fit
  res$pis <- pis
  res$taus <- taus
  return(res)
}

# replicate weights for non-ippm distributions -
"initiate_fit_sam" <- function(y, X, site_spp_weights, offset, y_is_na, G, S, disty, control){
  fm_sp_mods <- surveillance::plapply(1:S, ecomix:::apply_glmnet_sam, y, X, site_spp_weights, offset, y_is_na, disty,
                                   .parallel = control$cores, .verbose = !control$quiet)

  alphas <- unlist(lapply(fm_sp_mods, `[[`, 1))
  betas <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))
  disp <- unlist(lapply(fm_sp_mods, `[[`, 3))

  if(control$init_method=='kmeans'){
    if(!control$quiet)message( "Initial groups by K-means clustering\n")
    tmp1 <- stats::kmeans(betas, centers=G, nstart=100)
    tmp_grp <- tmp1$cluster
    grp_coefs <- apply(betas, 2, function(x) tapply(x, tmp_grp, mean))
  }

  if(control$init_method=='random' | is.null(tmp_grp)){
    if(!control$quiet)message( "Initial groups by random allocation and means from random numbers\n")
    grp_coefs <- matrix( stats::rnorm(G*ncol(betas), sd=control$init_sd, mean=0), nrow=G, ncol=ncol(betas))
    tmp_grp <- sample(1:G, S, replace=TRUE)
  }

  colnames(grp_coefs) <- colnames(X[,-1])
  results <- list()
  results$grps <- tmp_grp
  results$alphas <- alphas
  results$betas <- grp_coefs
  results$disp <- disp

  return(results)
}

"get_logls_sam" <- function(first_fit, fits, spp_weights, G, S, disty){

  if(is.null(spp_weights))spp_weights <- rep(1,S) #for bayesian boostrap.

  #bernoulli
  if(disty==1){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    link <- stats::make.link(link = "logit")
    for(ss in 1:S){
      for(gg in 1:G){
        lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(lp),log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #poisson
  if(disty==2){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss],exp(lp),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #ippm
  if(disty==3){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      sp_idx<-!first_fit$y_is_na[,ss]
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[sp_idx,1] * fits$alphas[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$betas[gg,] + first_fit$offset[sp_idx]
        logl_sp[ss,gg] <- first_fit$y[sp_idx,ss] %*% lp - first_fit$weights[sp_idx,ss] %*% exp(lp)
      }
    }
  }
  #negative binomial
  if(disty==4){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss],mu=exp(lp),size=fits$disp[ss],log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #tweedie
  if(disty==5){
  stop('no tweedie')

  }
  # gaussian
  if(disty==6){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=lp,sd=fits$disp[ss],log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  return(logl_sp)
}

fitmix_EM_sam <- function(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){


  ite <- 1
  logl_old <- -99999999
  logl_new <- -88888888

  # get starting values
  starting_values <- get_initial_values_sam(y = y, X = X,
                                            spp_weights = spp_weights,
                                            site_spp_weights = site_spp_weights,
                                            offset = offset, y_is_na = y_is_na,
                                            G = G, S = S,
                                            disty=disty,
                                            control = control)

  # first e-step
  pis <- starting_values$pis
  fits <- starting_values$fits
  taus <- starting_values$taus
  first_fit <- starting_values$first_fit

  while(control$em_reltol(logl_new,logl_old) & ite <= control$em_steps){

  # Estimate pis from tau.
  pis <- colSums(taus)/S

  if (any(pis == 0)) {
    starting_values <- get_initial_values_sam(y = y, X = X,
                                              spp_weights = spp_weights,
                                              site_spp_weights = site_spp_weights,
                                              offset = offset, y_is_na = y_is_na,
                                              G = G, S = S,
                                              disty=disty,
                                              control = control)
    pis <- starting_values$pis
    fits <- starting_values$fits
    taus <- starting_values$taus
    first_fit <- starting_values$first_fit
    ite <- 1
  }

  # m-step
  tmp <- stats::nlminb(start=fits$betas, objective=incom_logl_mix_coefs, gradient=NULL, hessian=NULL,
                       eta=additive_logistic(pis,inv = TRUE)[-G],
                       first_fit = first_fit,
                       fits = fits,
                       spp_weights = spp_weights,
                       G=G, S=S,
                       disty = disty)
  fits$betas <- update_mix_coefs(fits$betas, tmp$par)

  fm_sp_int <- surveillance::plapply(1:S, apply_glmnet_sam, y, X, site_spp_weights, offset, y_is_na, disty, .parallel = control$cores, .verbose = !control$quiet) #check weights in this.
  alphas <- unlist(lapply(fm_sp_int, `[[`, 1))
  fits$alphas <- update_sp_coefs(fits$alphas,alphas)

  # e-step
  # get the log-likes and taus
  logls <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- get_taus(pis, logls, G, S)

  #update the likelihood
  logl_old <- logl_new
  logl_new <- get_incomplete_logl_sam(eta = additive_logistic(pis,inv = TRUE)[-G], first_fit, fits, spp_weights, G, S, disty)
  ite <- ite + 1
  }

  taus <- data.frame(taus)
  names(taus) <- paste("grp.", 1:G, sep = "")
  int_out <- fits$alphas
  fm_out <- fits$betas
  names(pis) <- paste("G", 1:G, sep = ".")
  eta <- additive_logistic(pis, TRUE)[-1]

  # estimate log-likelihood
  logl_new <- get_incomplete_logl_bernoulli_sp_function(eta, first_fit, fits, G, S)

  return(list(logl = logl_new, alpha = int_out, beta = fm_out,
  eta = eta, pis = pis, taus = round(taus,4)))

}

"incom_logl_mix_coefs" <- function(x, eta, first_fit, fits, spp_weights, G, S, disty){

  fits$betas <- matrix(x,nrow=nrow(fits$betas),ncol=ncol(fits$betas))
  tmp <- get_incomplete_logl_sam(eta, first_fit, fits, spp_weights, G, S, disty)
  return(-tmp)
}

# this appears to be working.
"get_incomplete_logl_sam" <- function(eta, first_fit, fits, spp_weights, G, S, disty){

  pis <- additive_logistic(eta)
  if(is.null(spp_weights))spp_weights <- rep(1,S) #for bayesian boostrap.

  #bernoulli
  if(disty==1){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    link <- stats::make.link(link = "logit")
    for(ss in 1:S){
      for(gg in 1:G){
        lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(lp),log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #poisson
  if(disty==2){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss],exp(lp),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #ippm
  if(disty==3){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      sp_idx<-!first_fit$y_is_na[,ss]
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[sp_idx,1] * fits$alphas[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$betas[gg,] + first_fit$offset[sp_idx]
        logl_sp[ss,gg] <- first_fit$y[sp_idx,ss] %*% lp - first_fit$weights[sp_idx,ss] %*% exp(lp)
      }
    }
  }
  #negative binomial
  if(disty==4){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss],mu=exp(lp),size=fits$disp[ss],log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #tweedie
  if(disty==5){
    stop('no tweedie')

  }
  # gaussian
  if(disty==6){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=lp,sd=fits$disp[ss],log=TRUE))
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



