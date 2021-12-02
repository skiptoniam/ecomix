
#'@export
"glm.nbinom" <- function(form, data, weights=NULL, offset=NULL, mustart=NULL, est_var=FALSE){
    X <- stats::model.matrix(form, data)
    t1 <- stats::model.frame(form, data)
    y <- stats::model.response(t1)
    if(is.null(offset)) offset <- stats::model.offset(t1)
    if(is.null(offset)) offset <- rep(0,length(y))
    if(is.null(weights)) weights <- stats::model.weights(t1)
    if(is.null(weights)) weights <- rep(1,length(y))

    fit <- glm.fit.nbinom(X, y, offset, weights, mustart, est_var)

    return(fit)

}

#'@export
"glm.fit.nbinom" <- function(x, y, offset = NULL, weights = NULL,
                             mustart = NULL, est_var = FALSE){
  X <- x
  if (is.null(offset))
    offset <- 0
  if (is.null(weights))
    weights <- rep(1, length(y))
  gradient <- rep(0, ncol(X) + 1)
  if (is.null(mustart)) {
    pars <- gradient + 1
  }
  else {
    pars <- mustart
  }
  fitted.values <- rep(0, length(y))
  logl <- .Call("Neg_Bin", pars, X, y, weights, offset, gradient,
                fitted.values, PACKAGE = "nbglm")
  vcov <- 0
  se <- rep(0, length(pars))
  if (est_var) {
    calc_deriv <- function(p) {
      gradient <- rep(0, length(pars))
      ll <- .Call("Neg_Bin_Gradient", p, X, y, weights,
                  offset, gradient, PACKAGE = "nbglm")
      return(gradient)
    }
    hes <- numDeriv::jacobian(calc_deriv,pars)
    # hes <- nd2(pars, calc.deriv)
    dim(hes) <- rep(length(pars), 2)
    vcov <- try(solve(hes))
    se <- try(sqrt(diag(vcov)))
    colnames(vcov) <- rownames(vcov) <- c("theta", colnames(X))
  }
  names(pars) <- names(se) <- names(gradient) <- c("theta", colnames(X))
  return(list(logl = logl, coef = pars[-1], theta = pars[1],
              se = se[-1], se.theta = se[1], fitted = fitted.values,
              gradient = gradient, vcov = vcov))
}
