library('mgcv')
set.seed(20)
n <- 100
S <- 20
beta <- 4
G <- 4
sig <- 2

alpha <- runif(S,-7,-4) #intercept
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
f3 <- function(x) 0 * x
theta <-  log( 1 + rgamma( n=S, shape=1, scale=0.75)) # dispersion

x0 <- runif(n)
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)

X <- cbind(x0,x1,x2,x3)

fb1 <- f0(x0) + f1(x1) + f2(x2)
fb2 <- f0(x0) + f2(x1)
fb3 <- f1(x1) + f2(x2)
fb4 <- f1(x1) + f3(x3)
scale <- 1

library(nlme)
## simulate truth
# set.seed(1);n<-400;sig<-2
x <- 0:(n-1)/(n-1)
# f <- 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
## produce scaled covariance matrix for AR1 errors...
V <- corMatrix(Initialize(corAR1(.6),data.frame(x=x)))
Cv <- chol(V)  # t(Cv)%*%Cv=V
## Simulate AR1 errors ...
e <- t(Cv)%*%rnorm(n,0,sig) # so cov(e) = V * sig^2
## Observe truth + AR1 errors
y1 <- f0(x) + e
y2 <- f1(x) + e
y3 <- f2(x) + e
y4 <- f3(x) + e

par(mfrow=c(2,2))
plot(y1)
lines(f0(x))
plot(y2)
lines(f1(x))
plot(y3)
lines(f2(x))
plot(y4)
lines(f3(x))

f <- cbind(fb1,fb2,fb3,fb4)
fitted <- matrix(0, dim(X)[1], S)
group <- rep(0, S)
offset <- rep(0,dim(X)[1])
link <- make.link("log")
for (ss in seq_len(S)) {
  gg <- ceiling(stats::runif(1) * G)
  eta_spp <- alpha[ss] # could addin species specific gammas here (rather than just an intercept)
  eta_mix <- f[,gg]
  eta <- eta_spp + eta_mix + offset
  fitted[, ss] <- link$linkinv(eta*scale)
  group[ss] <- gg
}

outcomes <- matrix(rnbinom(n * S, mu=as.numeric( fitted), size=1/rep(exp(theta), each=n)), nrow = n, ncol = S)
pis <- tapply(group, group, length)/S

Y <- outcomes
y_is_na <- is.na(Y)
X <- X
W <- matrix(1,nrow(X),1)
offy <- offset
site_spp_weights <- wts <- matrix(1,nrow(X),S)
spp_weights <- rep(1,S)
colnames(Y) <- paste0("spp",seq_len(S))
archetype_formula <- as.formula(paste0('cbind(',paste(paste0('spp',1:S),collapse = ','),")~s(x0)+s(x1)+s(x2)+s(x3)"))
species_formula <- as.formula(~1)


df <- as.data.frame(cbind(Y,W,X))
fm<- gam(spp1~s(x0)+s(x1)+s(x2)+s(x3),data = df,family = "nb")

# Check that the theta from gam match size in dnbinom
logLik(fm)
sum(dnbinom(Y[,1], mu = fitted(fm), size = fm$family$getTheta(TRUE), log = TRUE))
# Xp <- predict(fm,type="lpmatrix")
# Xp%*%vcov(fm)%*%t(Xp)

# colnames(Y) <- paste0("spp",seq_len(S))

# forms <- list('archetype_formula'=archetype_formula,'species_formula'=species_formula)
# dat <- data.frame(Y,W,X)
# mm <- makeMatrices(forms,dat)


getGamMM <- function(fm,forms){
  rhs1 <- strsplit(deparse(forms[[1]][[3]]), " \\+ ")[[1]]
  rhs2 <- strsplit(deparse(forms[[2]][[2]]), " \\+ ")[[1]]
  mmX <- model.matrix(fm,exclude=rhs2)
  mmX <- mmX[,colSums(mmX)!=0][,-1,drop=FALSE] #drop intercept
  if(rhs2!="1") mmW <- model.matrix(fm,exclude=rhs1)
  else mmW <- model.matrix(fm)[,1,drop=FALSE]
  return(list(gamW=mmW,gamX=mmX))
}

merge.formula <- function(forms, ss, disty, ...){

  # get character strings of the names for the responses - should only be in the sam form.
  lhs1 <- forms[[1]][[2]]

  # get character strings of the right hand sides
  rhs1 <- strsplit(deparse(forms[[1]][[3]]), " \\+ ")[[1]]
  rhs2 <- strsplit(deparse(forms[[2]][[2]]), " \\+ ")[[1]]

  # create the merged rhs and lhs in character string form
  rhs <- c(rhs1, rhs2)
  lhs <- gsub("\\(","",lhs1)[-1]

  # reformulate function
  if(disty == 3 ) out <- reformulate(rhs, paste0(lhs[[ss]],"/wts"))
  else out <- reformulate(rhs, lhs[[ss]])

  environment(out) <- parent.frame()

  return(out)
}

## try fitting a bunch of single species GAMs.
apply_gamsam_inits <- function(ss, forms, Y, W, X, y_is_na, offy, wts, disty){

  # which family to use?
  if(disty == 1)
    fam <- binomial() #gam
  if(disty == 2 | disty == 3 | disty == 5)
    fam <- poisson()
  if(disty == 4)
    fam <- mgcv:::nb()
  if(disty == 6)
    fam <- gaussian()

  ## a bit of f#&king around with formulas, this will probably break.
  tmpform <- merge.formula(forms,ss,disty)
  df <- data.frame(Y,W,X)[!y_is_na[,ss],]  # remove NA sites if needs - ippm specific.
  fm <- mgcv::gam(tmpform, data =df, weights = wts[!y_is_na[,ss],ss], offset = offy, family=fam)
  # print(head(model.matrix(fm)))
  # cat("\n\n")

  ##estimate the starting dispersion parameter.
  my_coefs <- fm$coefficients

  # species intercpets
  alpha <- fm$coefficients[1]
  names(alpha) <- paste0("alpha.",colnames(Y[,ss,drop=FALSE]))
  # mixture coefs
  ## this could break if similar named coefs.
  beta <- my_coefs[grep(paste0(colnames(X),collapse = "|"),names(my_coefs))]
  # species coefs apart from intercept
  if(ncol(W)>1) gamma <-  my_coefs[grep(paste0(colnames(W),collapse = "|"),names(my_coefs))] else gamma <- -99999
  if(disty%in%4) theta <- fm$family$getTheta(TRUE)
  return(list(alpha = alpha, beta = beta, gamma = gamma, theta = theta, mm_gam=getGamMM(fm,forms)))

}

allspmods <- lapply(seq_len(S),apply_gamsam_inits, forms, Y, W, X, y_is_na, offy, wts, disty)

## wrapper around single species GAMs - which then groups the coefs.
initiate_fit_sam_gam <- function(forms, y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){


  fm_sp_mods <-  surveillance::plapply(seq_len(S), apply_gamsam_inits, forms,
                                       Y, W, X, y_is_na, offy, wts, disty,
                                       .parallel = control$cores, .verbose = FALSE)

  alpha <- unlist(lapply(fm_sp_mods, `[[`, 1))
  if(ncol(X)==1){
    beta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))[,-1,drop=FALSE]
  } else {
    beta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))
  }
  if(ncol(X)==0) beta <- rep(-999999,G)
  if(ncol(W)>1){
    gamma <- do.call(rbind,lapply(fm_sp_mods, `[[`, 3))
  } else {
    gamma <- unlist(lapply(fm_sp_mods, `[[`, 3))
  }

  theta <- unlist(lapply(fm_sp_mods, `[[`, 4))

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

  if(control$init_method=='random2'){
    fmmvnorm <- stats::kmeans(beta, centers=G, nstart=100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G)

    random_coefs <- ecomix:::sam_random_inits(alpha, grp_coefs, gamma, theta, S, G, X, W, disty, mult=0.3)
    alpha <- random_coefs[[1]]
    grp_coefs <- random_coefs[[2]]
    gamma <- random_coefs[[3]]
    if(disty%in%c(4,6))theta <- random_coefs[[4]]
  }

  if(ncol(beta[,,drop=FALSE])==1){
    names(grp_coefs)[2] <- names(beta)[2]
  } else {
    colnames(grp_coefs) <- colnames(beta)
  }

  #get taus as starting values
  if(G==1){
    taus <- matrix(1,nrow=ncol(y), ncol = G)
  } else {
    taus <- matrix(0,nrow=ncol(y), ncol= G)
    for(j in seq_len(S))taus[j,fmmvnorm$cluster[j]] <- 1

  }

  taus <- taus/rowSums(taus)
  taus <- ecomix:::shrink_taus(taus,G=G)
  pis <- colMeans(taus)

  results <- list()
  results$grps <- tmp_grp
  results$alpha <- alpha
  results$beta <- grp_coefs
  results$gamma <- gamma
  results$theta <- theta
  results$taus <- taus
  results$pis <- pis
  results$gamW <- fm_sp_mods[[1]]$mm_gam$gamW
  results$gamX <- fm_sp_mods[[1]]$mm_gam$gamX

  return(results)
}

startfits_samgam <- initiate_fit_sam_gam(forms, y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)

## Another wrapper around single species fit.
get_initial_values_samgam <- function(forms, y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){

  # get intial model fits
  starting_values <- initiate_fit_sam_gam(forms, y, X, W, spp_weights, site_spp_weights,
                                      offset, y_is_na, G, S, disty, control)

  S <- ncol(y)

  fits <- list(alpha=starting_values$alpha,
               beta=starting_values$beta,
               gamma=starting_values$gamma,
               theta=starting_values$theta)
  first_fit <- list(y = y, x = X, W = W, gamX = starting_values$gamX,
                    gamW = starting_values$gamW,
                    spp_weights = spp_weights,
                    site_spp_weights = site_spp_weights, offset = offset,
                    y_is_na = y_is_na,
                    removed_species = starting_values$species_to_remove)

  res <- list()
  res$fits <- fits
  res$first_fit <- first_fit
  res$taus <- starting_values$taus
  res$pis <- colSums(starting_values$taus)/S
  return(res)
}

get_logls_samgam <- function(first_fit, fits, spp_weights, G, S, disty, get_fitted=TRUE){

  if(get_fitted) fitted_values <- array(0,dim=c(G,nrow(first_fit$y),S))

  #bernoulli
  if(disty==1){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    link <- stats::make.link(link = "logit")
    for(ss in 1:S){
      for(gg in 1:G){
        if(ncol(first_fit$x)==1) eta <- rep(fits$alpha[ss],nrow(first_fit$x)) + first_fit$offset
        else eta <- fits$alpha[ss] + as.matrix(first_fit$gamX) %*% fits$beta[gg,] + first_fit$offset
        if(ncol(first_fit$W)>1) eta <- eta + as.matrix(first_fit$W[,-1]) %*% fits$gamma[ss,]
        if(get_fitted) fitted_values[gg,,ss] <- link$linkinv(eta)
        logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(eta),log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

    #poisson
    if(disty==2){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        if(ncol(first_fit$x)==1) eta <- rep(fits$alpha[ss],nrow(first_fit$x)) + first_fit$offset
        else eta <- fits$alpha[ss] + as.matrix(first_fit$x) %*% fits$beta[gg,] + first_fit$offset
        if(ncol(first_fit$W)>1) eta <- eta + as.matrix(first_fit$W[,-1]) %*% fits$gamma[ss,]
        if(get_fitted) fitted_values[gg,,ss] <- exp(eta)
        logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss],exp(eta),log=TRUE))
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
        if(ncol(first_fit$x)==1) eta <- rep(fits$alpha[ss],nrow(first_fit$x[sp_idx,])) + first_fit$offset[sp_idx,]
        else eta <- fits$alpha[ss] + as.matrix(first_fit$x[sp_idx,]) %*% fits$beta[gg,] + first_fit$offset[sp_idx]
        if(ncol(first_fit$W)>1) eta <- eta + as.matrix(first_fit$W[sp_idx,-1]) %*% fits$gamma[ss,]
        if(get_fitted) fitted_values[gg,sp_idx,ss] <- exp(eta)
        logl_sp[ss,gg] <- (first_fit$y[sp_idx,ss] %*% eta - first_fit$site_spp_weights[sp_idx,ss] %*% exp(eta))
      }
    }
    }

    #negative binomial
    if(disty==4){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #eta is the same as log_lambda (linear predictor)
        etaMix <- as.matrix(first_fit$gamX) %*% fits$beta[gg,] + first_fit$offset
        if(ncol(first_fit$gamW)==1) etaSpp <- as.matrix(first_fit$gamW) %*% fits$alpha[ss]
        else etaSpp <- s.matrix(first_fit$gamW) %*% c(fits$alpha[ss],fits$gamma[ss,])
        eta <- etaMix + etaSpp
        if(get_fitted) fitted_values[gg,,ss] <- exp(eta)
        logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss],mu=exp(eta),size=fits$theta[ss],log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
    }

    # gaussian
    if(disty==6){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #eta is the same as log_lambda (linear predictor)
        if(ncol(first_fit$x)==1) eta <- rep(fits$alpha[ss],nrow(first_fit$x)) + first_fit$offset
        else eta <- fits$alpha[ss] + as.matrix(first_fit$x) %*% fits$beta[gg,] + first_fit$offset
        if(ncol(first_fit$W)>1) eta <- eta + as.matrix(first_fit$W[,-1]) %*% fits$gamma[ss,]
        if(get_fitted) fitted_values[gg,,ss] <- eta
        logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=eta,sd=exp(fits$theta[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
    }
    out.list <- list(logl_sp=logl_sp)
    if(get_fitted) out.list$fitted = fitted_values
    return(out.list)
}


get_incomplete_logl_samgam <- function(eta, first_fit, fits,
                                      spp_weights, G, S, disty){

  pis <- ecomix:::additive_logistic(eta)
  if(is.null(spp_weights))spp_weights <- rep(1,S) #for bayesian boostrap.

  #bernoulli
  if(disty==1){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    link <- stats::make.link(link = "logit")
    for(ss in 1:S){
      for(gg in 1:G){
        eta <- fits$alpha[ss] + as.matrix(first_fit$x) %*% fits$beta[gg,] + first_fit$offset
        if(ncol(first_fit$W)>1) eta <- eta + as.matrix(first_fit$W[,-1]) %*% fits$gamma[ss,]
        logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(eta),log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #poisson
  if(disty==2){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    link <- stats::make.link(link = "log")
    for(ss in 1:S){
      for(gg in 1:G){
        eta <- fits$alpha[ss] + as.matrix(first_fit$x) %*% fits$beta[gg,] + first_fit$offset
        if(ncol(first_fit$W)>1) eta <- eta + as.matrix(first_fit$W[,-1,drop=FALSE]) %*% fits$gamma[ss,]
        logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss], lambda = link$linkinv(eta),log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #ippm
  if(disty==3){
    link <- stats::make.link(link = "log")
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      sp_idx<-!first_fit$y_is_na[,ss]
      for(gg in 1:G){
        #eta is the same as log_lambda (linear predictor)
        eta <- fits$alpha[ss] + as.matrix(first_fit$x[sp_idx,]) %*% fits$beta[gg,] + first_fit$offset[sp_idx]
        if(ncol(first_fit$W)>1) eta <- eta + as.matrix(first_fit$W[sp_idx,-1,drop=FALSE]) %*% fits$gamma[ss,]
        logl_sp[ss,gg] <- first_fit$y[sp_idx,ss] %*% eta - first_fit$site_spp_weights[sp_idx,ss] %*% exp(eta)
      }
    }
  }

  #negative binomial
  if(disty==4){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #eta is the same as log_lambda (linear predictor)
        etaMix <- as.matrix(first_fit$gamX) %*% fits$beta[gg,] + first_fit$offset
        if(ncol(first_fit$gamW)==1) etaSpp <- as.matrix(first_fit$gamW) %*% fits$alpha[ss]
        else etaSpp <- s.matrix(first_fit$gamW) %*% c(fits$alpha[ss],fits$gamma[ss,])
        eta <- etaMix + etaSpp
        if(get_fitted) fitted_values[gg,,ss] <- exp(eta)
        logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss],mu=exp(eta),size=fits$theta[ss],log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  # gaussian
  if(disty==6){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #eta is the same as log_lambda (linear predictor)
        eta <- fits$alpha[ss] + as.matrix(first_fit$x) %*% fits$beta[gg,] + first_fit$offset
        if(ncol(first_fit$W)>1) eta <- eta + as.matrix(first_fit$W[,-1,drop=FALSE]) %*% fits$gamma[ss,]
        logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=eta,sd=exp(fits$theta[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  ak <- logl_sp + matrix(rep(log( pis), each=S), nrow=S, ncol=G)
  am <- apply( ak, 1, max)
  ak <- exp( ak-am)
  sppLogls <- am + log( rowSums( ak))
  # theta.range=c(0.001, 10); pen.parm=1.25
  # if(disty==4)  sppLogls <- sppLogls  + dbeta((exp(-fits$theta[ss])-theta.range[1]) / (theta.range[2]- theta.range[1]), pen.parm, pen.parm, log=TRUE)
  logl <- sum( sppLogls)
  return(logl)
}


apply_glm_spp_coefs_samgam <- function(ss, y, X, W, G, taus,
                                       site_spp_weights,
                                       offset, y_is_na, disty, fits){
  # which family to use?
  if(disty == 1)
    fam <- binomial() #gam
  if(disty == 2 | disty == 3 | disty == 5)
    fam <- poisson()
  if(disty == 4)
    fam <- mgcv:::nb()
  if(disty == 6)
    fam <- gaussian()

  ids_i <- !y_is_na[,ss]

  if (disty==3){
    outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  } else {
    outcomes <- as.numeric(y[ids_i,ss])
  }
  out1 <- kronecker(rep( 1, G), outcomes)
  X1 <- kronecker(rep( 1, G), gamX[ids_i,,drop=FALSE])
  W1 <- kronecker(rep( 1, G), gamW[ids_i,,drop=FALSE])
  wts1 <- kronecker(rep( 1, G),
                    as.numeric(site_spp_weights[ids_i,ss]))*rep(taus[ss,],
                                                                each=length(site_spp_weights[ids_i,ss]))
  offy1 <- kronecker(rep( 1, G), offset[ids_i])
  offy2 <- gamX[ids_i,] %*% t(fits$beta)
  offy2 <- as.numeric(offy2)
  offy <- offy1 + offy2

  if(disty %in% c(1,2,3,6)){
    tmpform <- as.formula( paste('out1','-1+W1+offset(offy)', sep='~'))
    ft_sp <- try(mgcv::gam(tmpform, weights=wts1, family=fam))
    if (class(ft_sp) %in% 'try-error'){
      # print(paste0(ss,"\n"))
      my_coefs <- rep(NA, ncol(W1))
    } else {
      my_coefs <- stats::coef(ft_sp)
      names(my_coefs) <- c(colnames(y)[ss],colnames(W[,-1,drop=FALSE]))
    }
  }
  if(disty %in% c(4)){
    tmpform <- as.formula( paste('out1','-1+W1+offset(offy)', sep='~'))
    ft_sp <- try(mgcv::gam(tmpform, data = as.data.frame(W1), weights=wts1, family=mgcv::negbin(theta=fits$theta[ss])))
    kount1 <- 1
    while( any(class( ft_sp) %in% 'try-error') & kount1 < 10){
      kount1 <- kount1 + 1
      theta <- 10 * exp(-fits$theta[ss])
      ft_sp <- try(mgcv::gam(tmpform, weights=wts1, family=mgcv::negbin(theta=theta)))
    }
    my_coefs <- ft_sp$coefficients
    names(my_coefs) <- c(colnames(y)[ss],colnames(W[,-1,drop=FALSE]))
  }
  return(list(alpha = my_coefs[1], gamma = my_coefs[-1]))
}

apply_glm_mix_coefs_samgam <- function(gg, y, gamX, gamW, site_spp_weights, offset,
                                       y_is_na, disty,
                                       taus, fits, mus){

  ### setup the data stucture for this model.
  Y_taus <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
  X_no_NA <- list()
  W_no_NA <- list()
  for (jj in 1:ncol(y)){
    X_no_NA[[jj]] <- gamX[!y_is_na[,jj],,drop=FALSE]
    W_no_NA[[jj]] <- gamW[!y_is_na[,jj],,drop=FALSE]
  }
  X_taus <- do.call(rbind, X_no_NA)
  W_taus <- do.call(rbind, W_no_NA)
  n_ys <- sapply(X_no_NA,nrow)
  wts_taus <- rep(taus[,gg,drop=FALSE],c(n_ys))
  site_weights <- as.matrix(as.matrix(unlist(as.data.frame(site_spp_weights[!y_is_na]))))
  wts <- wts_taus*site_weights

  ## negative binomial weights
  # mus <- logls_mus$fitted
  if(disty==4) wts <- wts/(1+rep(fits$theta,n_ys)*as.vector(mus[gg,,]))

  offy_mat <- replicate(ncol(y),offset)
  offy1 <- unlist(as.data.frame(offy_mat[!y_is_na]))

  if(ncol(gamW)>1) {
    offy2 <- W %*% t(cbind(fits$alpha,fits$gamma))
  } else {
    offy2 <- W %*% c(fits$alpha)
  }
  offy2 <- unlist(as.data.frame(offy2[!y_is_na]))
  offy <- as.numeric(offy1 + offy2)

  # which family to use?
  if( disty == 1)
    fam <- binomial()
  if( disty %in% c(2,3,4,5))
    fam <- poisson()
  if( disty == 6)
    fam <- gaussian()


  if (disty==3){
    Y_taus <- as.matrix(Y_taus/site_weights)
  } else {
    Y_taus <- as.matrix(Y_taus)
  }

  if(disty %in% c(1,2,3,4,5,6)){ #don't use for tweedie - try and fit negative_binomial using glm.fit.nbinom
    tmpform <- as.formula( paste('Y_taus','-1+X_taus+offset(offy)', sep='~'))
    ft_mix <- try(mgcv::gam(tmpform, weights=wts, family=fam))
    if (class(ft_mix) %in% 'try-error'){
      mix_coefs <- rep(NA, ncol(X_taus))
    } else {
      mix_coefs <- ft_mix$coefficients
    }
  }

  return(c(mix_coefs))
}

starting_values <- get_initial_values_samgam(forms=forms, y = y, X = X, W = W,
                                          spp_weights = spp_weights,
                                          site_spp_weights = site_spp_weights,
                                          offset = offset, y_is_na = y_is_na,
                                          G = G, S = S,
                                          disty = disty,
                                          control = control)

apply_optimise_spp_theta_samgam <- function(ss, first_fit, fits,
                                       G, disty, pis,
                                       theta.range = c(0.0001,100)){
  thet <- optimise(f = theta.logl_samgam, interval = theta.range, ss, first_fit,
                   fits, G, disty, pis, theta.range,
                   maximum = TRUE)$maximum
  return(thet)
}

theta.logl_samgam <- function(theta, ss, first_fit, fits, G,
                         disty, pis, theta.range){

  if(disty==4)link <- make.link('log')
  if(disty==6)link <- make.link('identity')

  if(ncol(first_fit$W)>1){
    eta_spp <- first_fit$gamW %*% c(fits$alpha[ss],fits$gamma[ss,])
  }else{
    eta_spp <- first_fit$gamW %*% c(fits$alpha[ss])
  }
  eta_mix <- first_fit$gamX %*% t(fits$beta)
  offy <- first_fit$offset
  y <- first_fit$y[,ss]
  logls <- rep( NA, G)
  for( gg in 1:G){
    eta <- eta_spp + eta_mix[,gg] + offy
    mus <- link$linkinv(eta)
    if(disty==4)logls[gg] <- sum( dnbinom(x = y, mu = mus, size = theta, log=TRUE))
    if(disty==6)logls[gg] <- sum( dnorm(x = y, mean = mus, sd = theta, log=TRUE))
  }
  ak <- logls + log(pis)
  am <- max(ak)
  ak <- exp( ak-am)
  sppLogls <- am + log( sum( ak))

  # pen.max <- theta.range[2]
  # pen.min <- theta.range[1]
  # shape1 <- shape2 <- 1.25
  # if(disty==4)  sppLogls <- sppLogls + dbeta( (theta-pen.min) / (pen.max-pen.min), shape1, shape2, log=TRUE)

  return(sppLogls)
}


## fit ECM sam gam
fitmix_ECM_sam_gam <- function(forms, y, X, W, spp_weights, site_spp_weights, offset,
                               y_is_na, G, S, disty, control){

  ite <- 1
  restart_ite <- 1
  logl_old <- -99999999
  logl_new <- -88888888

  # get starting values
  starting_values <- get_initial_values_samgam(forms = forms, y = y, X = X, W = W,
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
  logls_mus <- get_logls_samgam(first_fit, fits, spp_weights, G, S, disty)
  init_steps <- 3
  gamX <- first_fit$gamX
  gamW <- first_fit$gamW

  while(control$em_reltol(logl_new,logl_old) & ite <= control$em_steps){
    if(restart_ite>10){
      message('cannot find good starting values with initialisation\n
              and random starting values\n\n
              Please check the number of groups and coefs.')
      break
    }

    pis <- colSums(taus)/S

    if (any(pis < sqrt(.Machine$double.eps))) {
      if(restart_ite==1){
        cat('Pis have gone to zero - restarting with new initialisation \n')
        starting_values <- get_initial_values_samgam(forms = forms, y = y, X = X, W = W,
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
        fits$beta <- matrix(rnorm(G*(ncol(X))),G,(ncol(X)))
        if(ncol(W)>1)fits$gamma <- matrix(rnorm(G*(ncol(W[,-1]))),S,(ncol(W[,-1])))
        fits$alpha <- rnorm(S)
        if (disty%in%c(4,6)) fits$theta <- rep(0.05,S)
      }
      ite <- 1
      restart_ite <- restart_ite + 1
    }

    # m-step
    fm_spp_coefs <- surveillance::plapply(seq_len(S),
                                          apply_glm_spp_coefs_samgam,
                                          y, gamX, gamW, G, taus, site_spp_weights,
                                          offset, y_is_na, disty, fits,
                                          .parallel = control$cores,
                                          .verbose = FALSE)

    # get and update the intercepts.
    alpha <- unlist(lapply(fm_spp_coefs, `[[`, 1))
    fits$alpha <- ecomix:::update_coefs(fits$alpha,alpha)

    # get and update the gamma if included in the model.
    if(ncol(W)>1) {
      gamma <- do.call(rbind,lapply(fm_spp_coefs, `[[`, 2))
      fits$gamma <- update_coefs(fits$gamma,gamma)
    } else {
      fits$gamma <- -99999
    }

    ## update the betas
    # if(disty==4){
    #   fm_mix_coefs <- nlminb(start=fits$beta, objective=incomplete_negbin_logl, gradient=NULL,
    #                          hessian=NULL, pis=pis, first_fit=first_fit, fits=fits, G=G, S=S)
    #   fm_mix_coefs_mat <- matrix(fm_mix_coefs$par,G,ncol(X))
    # } else {
    fm_mix_coefs <- surveillance::plapply(seq_len(G),
                                            apply_glm_mix_coefs_samgam,
                                            y, gamX, gamW, site_spp_weights,
                                            offset, y_is_na, disty, taus,
                                            fits, logls_mus$fitted,
                                            .parallel = control$cores,
                                            .verbose = FALSE)

    fm_mix_coefs_mat <- do.call(rbind,fm_mix_coefs)
    # }
    fits$beta <- ecomix:::update_coefs(fits$beta, fm_mix_coefs_mat)

    ## need a function here that updates the dispersion parameter.
    if(ite >= init_steps){
      if(disty%in%c(4)){
        # fits$theta <- exp(-fits$theta)
        fm_theta <- surveillance::plapply(seq_len(S), apply_optimise_spp_theta_samgam,
                                          first_fit, fits,
                                          G, disty, pis,
                                          .parallel = control$cores,
                                          .verbose = FALSE)
        theta <- unlist(lapply(fm_theta, `[[`, 1))
        # theta <- log(1/theta)
        fits$theta <- ecomix:::update_coefs(fits$theta,theta,control$update_kappa[2])
      }

      if(disty%in%c(6)){
        fits$theta <- exp(fits$theta)
        fm_theta <- surveillance::plapply(seq_len(S), apply_optimise_spp_theta,
                                          first_fit, fits,
                                          G, disty, pis,
                                          .parallel = control$cores,
                                          .verbose = FALSE)
        theta <- unlist(lapply(fm_theta, `[[`, 1))
        theta <- log(theta)
        fits$theta <- update_coefs(fits$theta,theta,control$update_kappa[2])
      }
    }

    ## up to e-step.
    # e-step - get the log-likes and taus
    logls_mus <- get_logls_samgam(first_fit, fits, spp_weights, G, S, disty)
    taus <- ecomix:::get_taus(pis, logls_mus$logl_sp, G, S)

    #update the likelihood
    logl_old <- logl_new
    logl_new <- get_incomplete_logl_samgam(eta = ecomix:::additive_logistic(pis,inv = TRUE)[-G],
                                        first_fit, fits, spp_weights, G, S, disty)
    if(!control$quiet)cat("Iteration ",ite,"\n")
    if(!control$quiet)cat("Loglike: ", logl_new,"\n")
    if(!control$quiet)cat("Pis: ", pis,"\n")
    ite <- ite + 1
  }

  taus <- data.frame(taus)
  names(taus) <- paste("grp.", 1:G, sep = "")
  names(pis) <- paste("G", 1:G, sep = ".")
  eta <- ecomix:::additive_logistic(pis, TRUE)[-G]

  return(list(logl = logl_new, alpha = fits$alpha, beta = fits$beta,
              gamma = fits$gamma, theta = fits$theta,
              eta = eta, pis = pis, taus = round(taus,4),
              first_fit = first_fit))

}

library(ecomix)
em_samgam <- fitmix_ECM_sam_gam(forms, y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty,control= species_mix.control(em_steps = 10))

archetype = ecomix:::sam_internal_pred_groups(alpha = em_samgam$alpha,
                                     beta = em_samgam$beta,
                                     gamma = em_samgam$gamma,
                                     taus = em_samgam$taus, G = G, S = S, X = gamX, W = gamW,
                                     offset = em_samgam$first_fit$offset, family = "negative_binomial")

df <- as.data.frame(cbind(Y,W,X))
fm<- gam(spp1~s(x0)+s(x1)+s(x2)+s(x3),data = df,family = "nb")
newd <- data.frame(x0=(0:100)/100,x1=(0:100)/100,x2=(0:100)/100,x3=(0:100)/100)
Xp <- predict(fm,newd,type="lpmatrix")

archetype = ecomix:::sam_internal_pred_groups(alpha = em_samgam$alpha,
                                              beta = em_samgam$beta,
                                              gamma = em_samgam$gamma,
                                              taus = em_samgam$taus, G = G, S = S, X = Xp[,-1], W = Xp[,1,drop=FALSE],
                                              offset = rep(0,nrow(Xp)), family = "negative_binomial")

matplot(log(archetype),type = 'l')

## these model matricies could be used for a TMB version.
makeMatrices <- function(forms, dat) {

  # results list
  res <- list()
  a <- model.frame(update.formula(forms[["archetype_formula"]], . ~ 1), data = dat)
  Y <- model.extract(a, "response")
  res$Y <- Y
  # remove intercept from archetype formula and whack in a dummy y
  sam.form <- update.formula(forms[["archetype_formula"]], y ~ -1 + .)
  # add intercept if missing and a dummy y
  spp.form <- update.formula(forms[["species_formula"]],y ~ 1 +. )

  ## whack in a dummy response variable
  dat$y <- dat[,1]
  gam_sam <- gam(sam.form, data = dat, method = "REML", fit = TRUE)
  # G <- gam(sam.form, data = dat, method = "REML", fit = FALSE)
  # gam.fit(G=G)
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
  mm <- predict(gam_spp, newdata = sppdat, type = "lpmatrix")
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

mm <- makeMatrices(forms, dat)






