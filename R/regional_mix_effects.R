#' @export

"effect_data.regional_mix" <- function(focal.predictors, mod, ngrid = 50, ...){

  focal.predictors <- c("x1.1","x1.2","x1.3","x2.1","x2.2","x2.3")
  mod <- fm_regional_mix
  ngrid <- 50

  if (is.null(mod$titbits))
    stop("Model doesn't contain all information required for effectsPlotData.
         Please supply model with titbits (from titbits=TRUE in regional_mix call)")

  Mode <- function(x, na.rm = FALSE) {
    if (na.rm) {
      x = x[!is.na(x)]
    }
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }

  ## set up the data objects
  X <- mod$titbits$X
  W <- mod$titbits$W

  ## set up the variables in the formula
  tt <- terms(mod$titbits$rcp_formula)
  tt <- delete.response(tt)
  vars <- all.vars(parse(text=tt))
  if(ncol(W)>1){
    tt2 <- terms(mod$titbits$species_formula)
    tt2 <- delete.response(tt2)
    vars <- c(vars,all.vars(parse(text=tt2)))
  }
  nvars = length(vars)

  pred.data <- mod$titbits$data
  pred.data <- pred.data[,vars]

  ## check for factors
  factors <- NULL
  for(ii in 1:nvars){
    factors[ii] <- is.factor(pred.data[,vars[ii]]) | is.character(pred.data[,vars[ii]])
  }

  ## check for binary data (i.e factors that have been expanded before modelling as per model.matrix)
  binary <- apply(pred.data,2,function(x) { all(x %in% 0:1) })

  ## check for focal.predictors in pred.data
  focal.ids <- lapply(focal.predictors, grep, colnames(pred.data))

  # lists for data structures
  mfs <- list() #catch model.frames
  f.focal <- list()
  v.focal <- list()
  n.focal <- list()
  for(i in 1:length(focal.ids)){

    f.focal[[i]] <- factors[focal.ids[[i]]]
    v.focal[[i]] <- pred.data[, focal.ids[[i]],drop=FALSE]
    n.focal[[i]] <- seq_len(nvars)[-unlist(focal.ids[[i]])]

    xx <- list()
    for(j in 1:length(focal.ids[[i]])){
      if(f.focal[[i]][j]) {
        xx[[j]] = levels(v.focal[[i]][j])
        ngrid = length(xx)
      } else {
        mi = min(v.focal[[i]][j])
        ma = max(v.focal[[i]][j])
        xx[[j]] = seq(mi, ma, length.out = ngrid)
      }
    }
    XDataNew = data.frame(xx, stringsAsFactors = TRUE)
    colnames(XDataNew) = vars[focal.ids[[i]]]
    for (k in seq_len(length(n.focal[[i]]))) {
      non.focal = n.focal[[i]][k]
      f.non.focal = factors[non.focal]
      v.non.focal = pred.data[, vars[non.focal]]
      if (f.non.focal) {
        XDataNew[, vars[non.focal]] = Mode(v.non.focal)
      }
      if (!f.non.focal) {
        v.non.focal = pred.data[, vars[non.focal]]
        XDataNew[, vars[non.focal]] = mean(v.non.focal)
      }
    }
    mfs[[i]] <- XDataNew[,vars]
  }

  names(mfs) <- focal.predictors
  class(mfs) <- "regional_mix_effectPlotData"

  return(mfs)

}
