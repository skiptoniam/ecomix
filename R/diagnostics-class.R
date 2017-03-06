#' @title diagnostics objects
#' @rdname diagnostics
#' @name AIC
#' @description diagnostic tools for checking \code{ecomix} models; see \code{ecomix} for more specifics about how
#' run an ecomix model.
#' @param object an ecomix object returned from \code{ecomix}.
#' @export

AIC.ecomix <-
  function (object, ..., k = 2)
  {
    p <- length(unlist(object$coefs))
    if (is.null(k))
      k <- 2
    star.ic <- -2 * object$logl + k * p
    return(star.ic)
  }


#' @name BIC
#' @export
BIC.ecomix <-
  function (object, ...)
  {
    p <- length(unlist(object$coefs))
    k <- log(object$n)
    star.ic <- -2 * object$logl + k * p
    return(star.ic)
  }
