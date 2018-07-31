## generic functions for species mix.

"lambda_penalisation_fun" <- function(x,lambda,kappa=0.1){ #assumes that x spans to pretty-well the unpenalised estiamtes
  min.effective.penalty <- min( which( abs( x-tail( x, 1)) < 0.01 * abs( tail( x, 1))))    #the first that lambda that gives a coef close to the last lambda's corresponding coef
  min.effective.penalty <- lambda[min.effective.penalty]
  target.penalty <- kappa * min.effective.penalty
  res.pos <- which.min( (lambda-target.penalty)^2)
  res <- x[res.pos]
  return( res)
}

"apply_glmnet_sam" <- function(ss, y, X, weights, offset, y_is_na, disty){

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
  disp <- NULL
  if( disty == 4){
    locat.s <- lambda.seq[max(which(as.matrix(glmnet::coef.glmnet(ft_sp))==my_coefs,arr.ind = TRUE)[,2])]
    preds <-as.numeric(predict(ft_sp, s=locat.s,
                               type="response",
                               newx=X[ids_i,-1],
                               offset=offset))
    tmp <- MASS::theta.mm(outcomes, preds,
                          weights=spp_weights,
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
                               offset=offset))
    disp <- log(sqrt(sum((outcomes - preds)^2)/length(outcomes)))  #should be something like the resid standard Deviation.
  }

  return(list(alpha = my_coefs[1], beta = my_coefs[-1], disp = disp))

}
