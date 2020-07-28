## firstly let's generate some data!
x <- y <- 1:100 / 100
grid2D <- expand.grid( x, y)
grid2D$cellArea <- rep( 1/200, nrow( grid2D))  #all cells have same size here
grid2D$x1 <- runif(nrow(grid2D))
grid2D$x2 <- runif(rnorm(grid2D))

n_sp <- 25
LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
sp_name <- LETTERS702[1:(n_sp)]
n_g <- 3

set.seed(123)
thetas <- as.matrix(data.frame(beta1=rnorm(n_g),beta2=rnorm(n_g)))
X <- as.matrix(data.frame(const=1,x1=grid2D$x1,x2=grid2D$x2))

lambdas <- matrix(0, dim(X)[1], n_sp, dimnames=list(NULL,sp_name))
sp_int <- rep(0, n_sp)
group <- rep(0, n_sp)
for (s in 1:n_sp) {
  g <- ceiling(runif(1) * n_g)
  sp_int[s] <- runif(1, -3, -1)
  log_lambda <-  X%*%c(sp_int[s],thetas[g,])
  lambdas[, s] <- exp(log_lambda)
  group[s] <- g
}

# cell_size <- rep(mean(res(preds)),nrow(lambdas))

LAMBDAS <- apply(lambdas,2,function(x)sum(x*grid2D$cellArea))
Ns <- sapply(LAMBDAS,function(x)rpois(n=1, lambda= x))  #the observed number of presences
preds_df <- data.frame(idx=1:nrow(X),X)
presences <- list()
for(i in seq_len(n_sp)){
  presences[[i]] <- sample(x=preds_df$idx,size=Ns[i], replace=FALSE, prob=lambdas[,i]/LAMBDAS[i])
}

presence_coords <- lapply(presences,function(x)grid2D[x,1:2])
presences_sort <- lapply(presences,sort)

sp_dat_po_ul<-data.frame(sp=rep(sp_name,unlist(lapply(presences,length))),cell_num=unlist(presences_sort))
po_matrix <- bbgdm::table2pam(sp_dat_po_ul,site.id='cell_num',sp.id='sp')
po_matrix[po_matrix==0]<-NA
po_covariates <- X[as.numeric(rownames(po_matrix)),]
presence_data <- data.frame(po_matrix,po_covariates)
bkdata <- cbind(matrix(0,nrow(X),n_sp),X)
colnames(bkdata) <- colnames(presence_data)
mm <- rbind(presence_data,bkdata)
y <- mm[,1:n_sp]
y_is_na <- is.na(y)
X <- mm[,(n_sp+1):(1+n_sp+2)]
X[,-1] <- scale(X[,-1])
wts <- as.matrix(y)#absence_weights for now let's make these all one.
wts[!is.na(wts)]<-1
offy <- rep(0,nrow(y))
G <- 3
S <- 25

library(ecomix)
sam_form <- as.formula(paste0('cbind(',paste(LETTERS702[1:(n_sp)],collapse = ','),")~1+x1+x2"))
sp_form <- ~ 1
model_data <- make_mixture_data(species_data = as.matrix(y),
                                covariate_data = as.matrix(X[,-1]))
fm1 <- species_mix(sam_form, sp_form, data = model_data,  weights = wts, distribution = 'ippm', n_mixtures=3)


## now let's generate the model matrix as per ippm model structure.
mm <- rbind(presence_data,bkdata)
y <- mm[,1:n_sp]
y_is_na <- is.na(y)
X <- mm[,(n_sp+1):(1+n_sp+2)]
X[,-1] <- scale(X[,-1])
wts <- y#absence_weights for now let's make these all one.
wts[]<-1
offy <- rep(0,nrow(y))
G <- 3
S <- 25
# control <- species_mix.control(nlminb_trace = 1)

# fitmix_ippm(y,X,G,weights,offset,y_is_na)

set.seed(213)
source('./developing_functions/species_mix_ppm_functions_version3.R')
starting_values <- get_initial_values_ippm(y, X, wts, offy, y_is_na, G, S, cores=1, inits='kmeans', init.sd=1)
pis <- starting_values$pis
fits <- starting_values$fits
taus <- starting_values$taus
first_fit <- starting_values$first_fit

incom_logl_ippm_alpha <- function(x, first_fit, tpis, fits, G, S){
  fits$sp_intercepts <- x
  pis <- additive_logistic(tpis)
  tmp <- get_incomplete_logl_ippm_function(pis, first_fit, fits, G, S)
  return(-tmp)
}

incom_logl_ippm_beta <- function(x, first_fit, tpis, fits, G, S){

  fits$mix_coefs <- matrix(x,nrow=nrow(fits$mix_coef),ncol=ncol(fits$mix_coef))
  pis <- additive_logistic(tpis)
  tmp <- get_incomplete_logl_ippm_function(pis, first_fit, fits, G, S)
  return(-tmp)
}

incom_logl_ippm_pi <- function(x, first_fit, fits, G, S){
  tpis <- x
  pis <- additive_logistic(tpis)
  tmp <- get_incomplete_logl_ippm_function(pis, first_fit, fits, G, S)
  return(-tmp)
}

get_incomplete_logl_ippm <-  function(tpis, first_fit, fits, G, S){

  pis <- additive_logistic(tpis)
  logl_sp_ippm <- matrix(NA,nrow=S, ncol=G)
  for(ss in 1:S){
    sp_idx<-!first_fit$y_is_na[,ss]
    for(gg in 1:G){
      #eta is the same as log_mu (linear predictor)
      eta <- first_fit$x[sp_idx,1] * fits$sp_intercepts[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$mix_coefs[gg,] + first_fit$offset[sp_idx]

      logl_sp_ippm[ss,gg] <- first_fit$y[sp_idx,ss] %*% eta - first_fit$weights[sp_idx,ss] %*% exp(eta)
    }
  }
  ak <- logl_sp_ippm + matrix(rep(log(pis), each=S), nrow=S, ncol=G)
  am <- apply( ak, 1, max)
  ak <- exp( ak-am)
  sppLogls <- am + log( rowSums( ak))
  logl <- sum( sppLogls)
  return( logl)
}

# Notes on how to calculate dfdalphaa
# sp_idx <- !y_is_na[,s]
# # df/dmu
# # f(y)=(e^(-mu)*mu^y)/y!
# # take the log
# # logf(y) = -log(y!)+ylog(mu)-mu
# # df/dmu = (y/mu)-1
# dlogfdm <- exp(y[sp_idx,s]/mu)-1)
# dfdm <- y[sp_idx,s]/mu)-1
#
# # dmu/deta
# # eta_i <- alpha_i + X%*%beta_k
# # mu <- exp(eta)
# dmde <- exp(eta)
#
# # deta/dalpha = alpha_j + XB_k
# # deta/dalpha = 1 #of n length.
# deda <- rep(1,n)
get_dlogfdalpha <-  function(first_fit, fits, G, S){
  dlogfdalpha <- matrix(NA,nrow=S, ncol=G)
  for(ss in 1:S){
    sp_idx<-!first_fit$y_is_na[,ss]
    for(gg in 1:G){
      eta <- first_fit$x[sp_idx,1] * fits$sp_intercepts[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$mix_coefs[gg,] + first_fit$offset[sp_idx]
      dlogfdalpha[ss,gg] <- c(c(((first_fit$y[sp_idx,ss]/first_fit$weights[sp_idx,ss])/exp(eta))-1) %*% exp(eta) * 1) ## this should be exp(eta) or mu, which are the same thing!
      #Replace y for z.
    }
  }
  return(dlogfdalpha)
}

get_logls_ippm<-function(first_fit, fits, G, S){
  logl_sp_ippm <- matrix(NA, nrow=S, ncol=G)
  for(ss in 1:S){
    sp_idx<-!first_fit$y_is_na[,ss]
    for(gg in 1:G){
      #eta is the same as log_lambda (linear predictor)
      eta <- first_fit$x[sp_idx,1] * fits$sp_intercepts[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$mix_coefs[gg,] + first_fit$offset[sp_idx]
      logl_sp_ippm[ss,gg] <- first_fit$y[sp_idx,ss] %*% eta - first_fit$weights[sp_idx,ss] %*% exp(eta)
      #Need weights*z*lp + weights*lp.  Should be OK?  Check!  SDF
    }
  }
  return(logl_sp_ippm)
}

get_spp_logls<-function(first_fit, fits, G, S){
  logl_sp_ippm <- matrix(NA, nrow=S, ncol=G)
  for(ss in 1:S){
    sp_idx<-!first_fit$y_is_na[,ss]
    for(gg in 1:G){
      #eta is the same as log_lambda (linear predictor)
      eta <- first_fit$x[sp_idx,1] * fits$sp_intercepts[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$mix_coefs[gg,] + first_fit$offset[sp_idx]
      cat(fits$sp_intercepts[ss],fits$mix_coefs[gg,],"\n\n")
      logl_sp_ippm[ss,gg] <- first_fit$y[sp_idx,ss] %*% eta - first_fit$weights[sp_idx,ss] %*% exp(eta)
    }
  }
  ak <- logl_sp_ippm + matrix(rep(log(pis), each=S), nrow=S, ncol=G)
  am <- apply( ak, 1, max)
  ak <- exp( ak-am)
  sppLogls <- am + log( rowSums( ak))
  return(sppLogls)
}

species_logl_each_pi <- get_logls_ippm(first_fit, fits, G, S)
spplogls <- get_spp_logls(first_fit, fits, G, S)
dlogfdalpha <- get_dlogfdalpha(first_fit, fits, G, S)

gradient_fun_alpha <- function(tpis, species_logl_each_pi, spplogls, dlogfdalpha, G, S){
  grads <- matrix(NA, nrow=S, ncol=G)
  pis <- additive_logistic(tpis)
  for(gg in 1:G){ #GO through all pi's
        for(ss in 1:S){
        #calculate for alphas
        grads[ss, gg] <- exp(-(spplogls[ss]) + species_logl_each_pi[ss,gg] + log(pis[gg])) * dlogfdalpha[ss,gg]#* (sum(ippm_weights[,ss]/species_logl_each_pi[ss,gg])
        }
  }
  -apply(grads,1,sum)# need -ve because of the way the loglike optim is written out.
}

## now let's try and work out the dbeta
get_dlogfdbeta <-  function(first_fit, fits, G, S){
  Nc <- ncol(first_fit$x[,-1])
  X <- first_fit$x[,-1]
  dlogfdbeta <- array(NA,dim = c(G,Nc,S))
   for(gg in 1:G){
       for(ii in 1:Nc){
           for(ss in 1:S){
               sp_idx<-!first_fit$y_is_na[,ss]
               eta <- first_fit$x[sp_idx,1] * fits$sp_intercepts[ss] + as.matrix(X[sp_idx,]) %*% fits$mix_coefs[gg,] + first_fit$offset[sp_idx]
               # c(c(((first_fit$y[sp_idx,ss]/first_fit$weights[sp_idx,ss])/exp(eta))-1) %*% exp(eta)
               dlogfdbeta[gg,ii,ss] <- c(c(drop(((first_fit$y[sp_idx,ss]/first_fit$weights[sp_idx,ss])/exp(eta))-1) * exp(eta)) %*% as.matrix(X[sp_idx,ii,drop=FALSE])) ## this should be exp(eta) or mu, which are the same thing!
      }
    }
  }
  return(dlogfdbeta)
}

species_logl_each_pi <- get_logls_ippm(first_fit, fits, G, S)
spplogls <- get_spp_logls(first_fit, fits, G, S)
dlogfdbeta <- get_dlogfdbeta(first_fit, fits, G, S)

gradient_fun_betas <- function(pis, species_logl_each_pi, spplogls, dlogfdbeta, G, S){
  Nc <- dim(dlogfdbeta)[2]
  grads_beta <- matrix(0, nrow=G,ncol=Nc)
  for(gg in 1:G){ #GO through all pi's
    for(ii in 1:Nc){
      for(ss in 1:S){
      #calculate for alphas
      grads_beta[gg,ii] <- grads_beta[gg,ii] + exp(-(spplogls[ss]) + (species_logl_each_pi[ss,gg]) + log(pis[gg])) * dlogfdbeta[gg,ii,ss]
      }
    }
  }
  -c(grads_beta)# need -ve because of the way the loglike optim is written out.
}

## let's try and calculate the pi derivates
## this should calculate the derivate w.r.t pi
gradient_fun_pi <- function(tpis, species_logl_each_pi, spplogls, G, S){

  npis <- additive_logistic(tpis)
  dlogdpi <- matrix(0, nrow=G,ncol=S)
  for(gg in 1:G){ #GO through all pi's
      for(ss in 1:S){
        #calculate dlogdpi
        dlogdpi[gg,ss] <- exp(-(spplogls[ss]) + (species_logl_each_pi[ss,gg]))
      }
    }
  dlogdpis <- rowSums(dlogdpi)

  add_log_trans <-0
  pi_mat_deriv <- matrix(0,(G-1),(G))
  piDerivs <- rep(0,G-1)

  # for(g in seq_len(G-1)){
  #   add_log_trans = add_log_trans + exp(tpis[g])
  # }

  add_log_trans <- sum(exp(tpis)) + 1

  # add_log_trans = 1 + add_log_trans

  ## fix indexing for pi_mat_deriv
  for(g in seq_len(G)){
    for(i in seq_len(G-1)){ #// go through eta's
          if(g<(G)){
            if(i==g){
            pi_mat_deriv[i,g] = (exp(tpis[i])/add_log_trans) - (exp(2*tpis[i])/(add_log_trans*add_log_trans))#;// diag
            pi_mat_deriv[i,G] = pi_mat_deriv[i,G] + pi_mat_deriv[i,g]#pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1)));
            }else{
            pi_mat_deriv[i,g] = (-1*(exp(tpis[i])*exp(tpis[g]))) / (add_log_trans*add_log_trans)#; //off-diag
            pi_mat_deriv[i,G] = pi_mat_deriv[i,G] + pi_mat_deriv[i,g]#;
        }
          }
      print(pi_mat_deriv[i,g])
    }
  }

  for(i in seq_len(G-1)) pi_mat_deriv[i,G] = pi_mat_deriv[i,G] * -1

  # pi_mat_deriv = pi_mat_deriv * -1

  for(i in seq_len(G-1)){
      for(g in seq_len(G)){
          piDerivs[i] = piDerivs[i] + (dlogdpis[g] * pi_mat_deriv[i,g])
          }
  }
  return(-piDerivs)
}


library(numDeriv)
#this should return the gradient for the alpha and betas.
tpis <- additive_logistic(pis,TRUE)[-G]
dfdalpha_num <- grad(incom_logl_ippm_alpha, fits$sp_intercepts, first_fit=first_fit, tpis=tpis, fits=fits, G=G, S=S)
dfdalpha_anl <- gradient_fun_alpha(tpis, species_logl_each_pi, spplogls, dlogfdalpha, G, S)
all.equal(dfdalpha_anl,dfdalpha_num)

dfdbeta_num <- grad(incom_logl_ippm_beta, fits$mix_coefs, first_fit=first_fit, pis=pis, fits=fits, G=G, S=S)
dfdbeta_anl  <- gradient_fun_betas(pis, species_logl_each_pi, spplogls, dlogfdbeta, G, S)
all.equal(dfdbeta_num,dfdbeta_anl)

dfdpi_num <- grad(incom_logl_ippm_pi, tpis, first_fit=first_fit, fits=fits, G=G, S=S)
dfdpi_anl <- gradient_fun_pi(tpis, species_logl_each_pi, spplogls, G, S)
all.equal(dfdpi_num,dfdpi_anl)

