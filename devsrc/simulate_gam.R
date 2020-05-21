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
y1 <- f0(x) + f1(x) + f2(x) + f3(x) + e
y2 <- f1(x) + e
y3 <- f2(x) + e
y4 <- f3(x) + e

par(mfrow=c(2,2))
plot(y1)
lines(f0(x) + f1(x) + f2(x) + f3(x))
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
  eta_spp <- alpha[ss]
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
sites_spp_weights <- matrix(1,nrow(X),S)
spp_wts <- rep(1,S)

df <- as.data.frame(cbind(Y,W,X))
# fm<- gam(y~s(x0)+s(x1),data = df,family = "nb")

colnames(Y) <- paste0("spp",seq_len(S))
archetype_formula <- as.formula(paste0('cbind(',paste(paste0('spp',1:S),collapse = ','),")~s(x0)+s(x1)+s(x2)+s(x3)"))
species_formula <- as.formula(~1)

forms <- list('archetype_formula'=archetype_formula,'species_formula'=species_formula)
dat <- data.frame(Y,W,X)
mm <- makeMatrices(forms,dat)

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


apply_gamsam_inits <- function(ss, forms, Y, W, X, y_is_na, offy, wts, disty){

  # which family to use?
  if(disty == 1)
    fam <- binomial() #gam
  if(disty == 2 | disty == 3 | disty == 5)
    fam <- poisson()
  if(disty == 4)
    fam <- nb()
  if(disty == 6)
    fam <- gaussian()

  ## a bit of f#&king around with formulas, this will probably break.
  tmpform <- merge.formula(forms,ss,disty)
  df <- data.frame(Y,W,X)[!y_is_na[,ss],]  # remove NA sites if needs - ippm specific.
  fm <- mgcv::gam(tmpform, data =df, weights = wts, offset = offy, family=fam)

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
  return(list(alpha = alpha, beta = beta, gamma = gamma, theta = theta))

}

## this appears to work :)
allspmods <- lapply(seq_len(S), initial_values_gamsam, forms, Y, W, X, y_is_na, offy, wts, disty)

initiate_fit_sam_gam <- function(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){


  fm_sp_mods <-  surveillance::plapply(seq_len(S), apply_gamsam_inits, forms,
                                       Y, W, X, y_is_na, offy, wts, disty,
                                       .parallel = control$cores, .verbose = FALSE)
  fm_sp_mods <-  surveillance::plapply(seq_len(S), apply_glmnet_sam_inits, y, X, W,
                                       site_spp_weights, offset, y_is_na, disty,
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

  species_to_remove <- which(apply(beta, 1, function(x) all(is.na(x))))

  if(length(species_to_remove)>0){
    #update fits
    alpha <- alpha[-species_to_remove]
    beta <- beta[-species_to_remove,,drop=FALSE]
    theta <- theta[-species_to_remove]

    # update y, y_is_na and weights
    updated_y <- update_species_data_structure(y, y_is_na, spp_weights, site_spp_weights, species_to_remove)
    y <- updated_y[[1]]
    y_is_na <- updated_y[[2]]
    site_spp_weights <- updated_y[[3]]
  } else {
    species_to_remove <- NA
  }

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

  if(control$init_method=='random2'){
    fmmvnorm <- stats::kmeans(beta, centers=G, nstart=100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G)

    random_coefs <- sam_random_inits(alpha, grp_coefs, gamma, theta, S, G, X, W, disty, mult=0.3)
    alpha <- random_coefs[[1]]
    grp_coefs <- random_coefs[[2]]
    gamma <- random_coefs[[3]]
    if(disty%in%c(4,6))theta <- random_coefs[[4]]
  }

  if(ncol(X[,,drop=FALSE])==1){ names(grp_coefs)[2] <- names(X[,,drop=FALSE])[2]
  } else {
    colnames(grp_coefs) <- colnames(X[,,drop=FALSE])
  }

  #get taus as starting values
  if(G==1){
    taus <- matrix(1,nrow=ncol(y), ncol = G)
  } else {
    taus <- matrix(0,nrow=ncol(y), ncol= G)
    if(length(sel_omit_spp)>0){
      for(j in 1:length((1:S)[-sel_omit_spp]))
        taus[(1:S)[-sel_omit_spp][j],fmmvnorm$cluster[j]] <- 1
      taus[sel_omit_spp,] <- matrix(runif(length(sel_omit_spp)*G),length(sel_omit_spp), G)
    } else {
      for(j in seq_len(S))taus[j,fmmvnorm$cluster[j]] <- 1
    }
  }

  taus <- taus/rowSums(taus)
  taus <- shrink_taus(taus,G=G)
  pis <- colMeans(taus)

  results <- list()
  results$grps <- tmp_grp
  results$alpha <- alpha
  results$beta <- grp_coefs
  results$gamma <- gamma
  results$theta <- theta
  results$taus <- taus
  results$pis <- pis
  results$species_to_remove <- species_to_remove

  return(results)
}


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

y <- mm$Y
X <- mm$A_sam

# y, X, W
# replace with mm from gam
disty <- 4

apply_gam_sam_inits <- function(ss, mm, site_spp_weights, offset, y_is_na, disty){

  # mm - this is the model matrix with all the include data.
  y <- mm$Y

  if(is.null(mm$Xs_sam)){
    X <- mm$A_sam
  } else {
    X <- cbind(mm$Xs_sam,mm$A_sam)
  }

  if(is.null(mm$Xs_spp)){
    W <- mm$A_spp
  } else {
    W <- cbind(mm$Xs_spp,mm$A_spp)
  }

  # which family to use?
  if(disty == 1)
    fam <- binomial() #gam
  if(disty == 2 | disty == 3 | disty == 5)
    fam <- poisson()
  if(disty == 4)
    fam <- nb()
  if(disty == 6)
    fam <- gaussian()

  ids_i <- !y_is_na[,ss]

  if (disty==3){
    outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  } else {
    outcomes <- as.matrix(y[ids_i,ss])
  }
#
#   if(ncol(X)==1){
#     X<-cbind(1,X[ids_i,,drop=FALSE])
#   }

  df <- cbind(W,X)


  if( disty %in% c(1,2,3,4,6)){
    # lambda.seq <- sort( unique( c( seq( from=1/0.001, to=1, length=25), seq( from=1/0.1, to=1, length=10))), decreasing=TRUE)

    ft_sp <- try(mgcv::gam.fit3(y=data.frame(outcomes), x=data.frame(df), sp = mm$gam_sam$sp, mm$S_sam, family = fam))
    if (any(class(ft_sp) %in% 'try-error')){
      my_coefs <- rep(NA, ncol(X[ids_i,]))
      names(my_coefs) <- colnames(cbind(X[ids_i,,drop=FALSE],W[ids_i,,drop=FALSE]))
    } else {
      if(ncol(X)==1) my_coefs <- t(as.matrix(my.coefs[-1]))
      my_coefs <- t(as.matrix(my.coefs))
    }

    ##estimate the starting dispersion parameter.
    theta <- -99999
    if( disty == 4){
      tmp <- MASS::theta.mm(outcomes, as.numeric(predict(ft_sp, s=locat.s,
                                                         type="response",
                                                         newx=as.matrix(df),
                                                         newoffset=offset[ids_i])),
                            weights=as.matrix(site_spp_weights[ids_i,ss]),
                            dfr=length(outcomes), eps=1e-4)
      if(tmp>2) tmp <- 2
      theta <- log( 1/tmp)
      # cat(tmp,"\n")
    }
    if( disty == 6){
      preds <- as.numeric( predict(ft_sp, s=locat.s, type="link",
                                   newx=as.matrix(df), newoffset=offset[ids_i]))
      theta <- log( sqrt( sum((outcomes - preds)^2)/length(outcomes)))  #should be something like the resid standard
    }
  }
  # species intercpets
  alpha <- my_coefs[1]
  # mixture coefs
  beta <- my_coefs[match(colnames(X), colnames(my_coefs))]
  # species coefs apart from intercept
  if(ncol(W)>1) gamma <-  my_coefs[match(colnames(W[,-1,drop=FALSE]), colnames(my_coefs))] else gamma <- -99999

  return(list(alpha = alpha, beta = beta, gamma = gamma, theta = theta))
}

initiate_fit_sam_gam <- function(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){


  fm_sp_mods <-  surveillance::plapply(seq_len(S), apply_glmnet_sam_inits, y, X, W,
                                       site_spp_weights, offset, y_is_na, disty,
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

  species_to_remove <- which(apply(beta, 1, function(x) all(is.na(x))))

  if(length(species_to_remove)>0){
    #update fits
    alpha <- alpha[-species_to_remove]
    beta <- beta[-species_to_remove,,drop=FALSE]
    theta <- theta[-species_to_remove]

    # update y, y_is_na and weights
    updated_y <- update_species_data_structure(y, y_is_na, spp_weights, site_spp_weights, species_to_remove)
    y <- updated_y[[1]]
    y_is_na <- updated_y[[2]]
    site_spp_weights <- updated_y[[3]]
  } else {
    species_to_remove <- NA
  }

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

  if(control$init_method=='random2'){
    fmmvnorm <- stats::kmeans(beta, centers=G, nstart=100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G)

    random_coefs <- sam_random_inits(alpha, grp_coefs, gamma, theta, S, G, X, W, disty, mult=0.3)
    alpha <- random_coefs[[1]]
    grp_coefs <- random_coefs[[2]]
    gamma <- random_coefs[[3]]
    if(disty%in%c(4,6))theta <- random_coefs[[4]]
  }

  if(ncol(X[,,drop=FALSE])==1){ names(grp_coefs)[2] <- names(X[,,drop=FALSE])[2]
  } else {
    colnames(grp_coefs) <- colnames(X[,,drop=FALSE])
  }

  #get taus as starting values
  if(G==1){
    taus <- matrix(1,nrow=ncol(y), ncol = G)
  } else {
    taus <- matrix(0,nrow=ncol(y), ncol= G)
    if(length(sel_omit_spp)>0){
      for(j in 1:length((1:S)[-sel_omit_spp]))
        taus[(1:S)[-sel_omit_spp][j],fmmvnorm$cluster[j]] <- 1
      taus[sel_omit_spp,] <- matrix(runif(length(sel_omit_spp)*G),length(sel_omit_spp), G)
    } else {
      for(j in seq_len(S))taus[j,fmmvnorm$cluster[j]] <- 1
    }
  }

  taus <- taus/rowSums(taus)
  taus <- shrink_taus(taus,G=G)
  pis <- colMeans(taus)

  results <- list()
  results$grps <- tmp_grp
  results$alpha <- alpha
  results$beta <- grp_coefs
  results$gamma <- gamma
  results$theta <- theta
  results$taus <- taus
  results$pis <- pis
  results$species_to_remove <- species_to_remove

  return(results)
}



"get_initial_values_sam" <- function(y, X, W, spp_weights, site_spp_weights,
                                     offset, y_is_na, G, S, disty, control){

  # get intial model fits
  starting_values <- initiate_fit_sam(y, X, W, spp_weights, site_spp_weights,
                                      offset, y_is_na, G, S, disty, control)

  #if any are errors then remove them from the models forever.
  # updated_y <- update_species_data_structure(y, y_is_na,
  #                                            spp_weights,
  #                                            site_spp_weights,
  #                                            starting_values$species_to_remove)
  # y <- updated_y[[1]]
  # y_is_na <- updated_y[[2]]
  # spp_weights <- updated_y[[3]]
  # site_spp_weights <- updated_y[[4]]
  S <- ncol(y)

  fits <- list(alpha=starting_values$alpha,
               beta=starting_values$beta,
               gamma=starting_values$gamma,
               theta=starting_values$theta)
  first_fit <- list(y = y, x = X, W = W, spp_weights = spp_weights,
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


fitmix_ECM_sam_gam <- function(y, X, W, spp_weights, site_spp_weights, offset,
                             y_is_na, G, S, disty, control){

  ite <- 1
  restart_ite <- 1
  logl_old <- -99999999
  logl_new <- -88888888

  # get starting values
  starting_values <- get_initial_values_sam(y = y, X = X, W = W,
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
  init_steps <- 3

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
        starting_values <- get_initial_values_sam(y = y, X = X, W = W,
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
                                          apply_glm_spp_coefs_sams,
                                          y, X, W, G, taus, site_spp_weights,
                                          offset, y_is_na, disty, fits,
                                          .parallel = control$cores,
                                          .verbose = FALSE)

    # get and update the intercepts.
    alpha <- unlist(lapply(fm_spp_coefs, `[[`, 1))
    fits$alpha <- update_coefs(fits$alpha,alpha)

    # get and update the gamma if included in the model.
    if(ncol(W)>1) {
      gamma <- do.call(rbind,lapply(fm_spp_coefs, `[[`, 2))
      fits$gamma <- update_coefs(fits$gamma,gamma)
    } else {
      fits$gamma <- -99999
    }

    ## update the betas
    if(disty==4){
      fm_mix_coefs <- nlminb(start=fits$beta, objective=incomplete_negbin_logl, gradient=NULL,
                             hessian=NULL, pis=pis, first_fit=first_fit, fits=fits, G=G, S=S)
      fm_mix_coefs_mat <- matrix(fm_mix_coefs$par,G,ncol(X))
    } else {
      fm_mix_coefs <- surveillance::plapply(seq_len(G),
                                            apply_glm_mix_coefs_sams,
                                            y, X, W, site_spp_weights,
                                            offset, y_is_na, disty, taus,
                                            fits, logls_mus$fitted,
                                            .parallel = control$cores,
                                            .verbose = FALSE)

      fm_mix_coefs_mat <- do.call(rbind,fm_mix_coefs)
    }
    fits$beta <- update_coefs(fits$beta, fm_mix_coefs_mat)

    ## need a function here that updates the dispersion parameter.
    if(ite >= init_steps){
      if(disty%in%c(4)){
        fits$theta <- exp(-fits$theta)
        fm_theta <- surveillance::plapply(seq_len(S), apply_optimise_spp_theta,
                                          first_fit, fits,
                                          G, disty, pis,
                                          .parallel = control$cores,
                                          .verbose = FALSE)
        theta <- unlist(lapply(fm_theta, `[[`, 1))
        theta <- log(1/theta)
        fits$theta <- update_coefs(log(1/fits$theta),theta,control$update_kappa[2])
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
        fits$theta <- update_coefs(log(fits$theta),theta,control$update_kappa[2])
      }
    }
    # e-step - get the log-likes and taus
    logls_mus <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
    taus <- get_taus(pis, logls_mus$logl_sp, G, S)

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
  names(pis) <- paste("G", 1:G, sep = ".")
  eta <- additive_logistic(pis, TRUE)[-G]

  return(list(logl = logl_new, alpha = fits$alpha, beta = fits$beta,
              gamma = fits$gamma, theta = fits$theta,
              eta = eta, pis = pis, taus = round(taus,4),
              first_fit = first_fit))

}




