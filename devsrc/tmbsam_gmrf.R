require( RandomFields )
location <- NULL
nsites <- 100
nspecies <- 50
SpatialScale=0.4; SD_A=0.5;
if(is.null(location)) locations <- cbind( "x"=runif(nsites, min=0,max=1), "y"=runif(nsites, min=0,max=1) )
model_A <- RandomFields::RMgauss(var=SD_A^2, scale=SpatialScale)
logMeanDens <- 1
alpha_p = (diag(nspecies)) %*% rep(logMeanDens,nspecies)

archetype_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
species_form <- as.formula(~1+w1+w2)
beta <- matrix(c(-3.6,0.5,
                 -0.9,1.0,
                 0.9,-2.9,
                 2.2,5.4),
               4,2,byrow=TRUE)
gamma <- matrix(c(rnorm(50,1),rnorm(50,-2)),50,2)
dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
dat[,-1] <- scale(dat[,-1])

#update the formula to old format.
if(!is.null(archetype_formula))
  archetype_formula <- stats::as.formula(archetype_formula)
if(!is.null(species_formula))
  species_formula <- stats::as.formula(species_formula)

sam_org <- archetype_formula
spp_org <- species_formula

# how many species to simulate???
S <- length(archetype_formula[[2]])-1
G <- n_mixtures

archetype_formula <- update(archetype_formula,y~.-1)
archetype_formula[[2]] <- NULL

X <- stats::model.matrix(archetype_formula, dat)
W <- stats::model.matrix(species_formula, dat)
if(is.null(offset)) offset <- rep(0,nrow(X))

n <- nrow(X)
npx <- ncol(X)
npw <- ncol(W) - 1

alpha <- NULL
if (is.null(alpha)) {
  message("Random alpha from normal (-1,0.5) distribution")
  alpha <- rnorm(S,-1,0.5)
}
if (is.null(beta) | length(beta) != (npx) * G) {
  message("Random values for beta")
  beta <- rnorm(npx*(G))
}
beta <- matrix(as.numeric(beta), nrow = G)
if(( is.null(gamma) | length( gamma) != S * npw)){
  if( npw != 0){
    message("Random values for gamma")
    gamma <- rnorm( npw*S)
    gamma <- matrix( as.numeric(gamma), nrow=S, ncol=npw)
  } else {
    gamma <- NULL
  }
} else {
  gamma <- matrix( as.numeric( gamma), nrow=S)
}
distribution <- "poisson"
theta <- NULL
if( distribution == "negative_binomial" & (is.null(theta) | length( theta) != S)){
  message( "Random values for overdispersions")
  theta <- log( 1 + rgamma( n=S, shape=1, scale=0.75))
}
if( distribution=="gaussian" & (is.null( theta) | length( theta) != S)){
  message( "Random values for species' variance parameters")
  theta <- log( 1 + rgamma( n+S, shape=1, scale=0.75))
}

if(distribution %in% 'bernoulli') link <- make.link('logit')
if(distribution %in% c('poisson','ippm','negative_binomial')) link <- make.link('log')
if(distribution %in% c('gaussian')) link <- make.link('identity')
if(distribution %in% 'ippm') {
  grid <- simulate_ippm_grid(X,W)
  grid2D <- grid$grid2D
  X <- grid$X
  W <- grid$W
}

logMeanDens <- 1
alpha_p = (diag(nspecies)) %*% rep(logMeanDens,nspecies)
A_sp = matrix(NA, nrow=nsites, ncol=nspecies)
for(p in 1:nspecies){
  A_sp[,p] = RandomFields::RFsimulate(model=model_A, x=locations[,'x'], y=locations[,'y'])@data[,1]
  A_sp[,p] = A_sp[,p] - mean(A_sp[,p]) + alpha_p[p]
}

logMeanDens <- 1
beta_p = (diag(G)) %*% rep(logMeanDens,G)
A_G = matrix(NA, nrow=nsites, ncol=G)
for(p in 1:G){
  A_G[,p] = RandomFields::RFsimulate(model=model_A, x=locations[,'x'], y=locations[,'y'])@data[,1]
  A_G[,p] = A_G[,p] - mean(A_G[,p]) + beta_p[p]
}


## simulate the groups and fitted values.
fitted <- matrix(0, dim(X)[1], S)
group <- rep(0, S)
for (ss in seq_len(S)) {
  gg <- ceiling(stats::runif(1) * G)
  eta_spp <- W %*% c(alpha[ss],gamma[ss,])
  eta_mix <- X %*% beta[gg, ]
  eta <- eta_spp + eta_mix +
  fitted[, ss] <- link$linkinv(eta)
  group[ss] <- gg
}

if( distribution=="bernoulli")
  outcomes <- matrix(rbinom(n * S, 1, as.numeric( fitted)), nrow = n, ncol = S)
if( distribution=="poisson")
  outcomes <- matrix(rpois(n * S, lambda=as.numeric( fitted)), nrow = n, ncol = S)
if( distribution=="ippm")
  outcomes <- simulate_ippm_outcomes(X, W, S, grid2D, fitted)
if( distribution=="negative_binomial")
  outcomes <- matrix(rnbinom(n * S, mu=as.numeric( fitted), size=1/rep(exp(theta), each=n)), nrow = n, ncol = S)
if( distribution=="gaussian")
  outcomes <- matrix( rnorm( n=n*S, mean=as.numeric( fitted), sd=rep( exp(theta), each=n)), nrow=n, ncol=S)

pi <- tapply(group, group, length)/S

if (distribution=='ippm'){
  res <- outcomes$mm
  wts <- outcomes$weights
  colnames(wts) <- all.vars(sam_org)[1:S]
  colnames(fitted) <- all.vars(sam_org)[1:S]
} else {
  colnames(fitted) <- all.vars(sam_org)[1:S]
  colnames(outcomes) <-  all.vars(sam_org)[1:S]
  if(ncol(W)>1){
    if( !all( offset==0))
      res <- data.frame(outcomes,const=1, X, W[,-1,drop=FALSE], offset=offset)
    else
      res <- data.frame(outcomes,const=1, X, W[,-1,drop=FALSE])
  } else {
    if( !all( offset==0))
      res <- data.frame(outcomes,const=1, X, offset=offset)
    else
      res <- data.frame(outcomes,const=1, X)
  }
  wts <- NULL
}
attr(res, "SAMs") <- group
attr(res, "pis") <- pi
attr(res, "alpha") <- alpha
attr(res, "beta") <- beta
attr(res, "gamma") <- gamma
attr(res, "theta") <- theta
attr(res, "mu") <- fitted
attr(res, "ippm_weights") <- wts


A_sp = matrix(NA, nrow=nsites, ncol=nspecies)
for(p in 1:nspecies){
  A_sp[,p] = RandomFields::RFsimulate(model=model_A, x=locations[,'x'], y=locations[,'y'])@data[,1]
  A_sp[,p] = A_sp[,p] - mean(A_sp[,p]) + alpha_p[p]
}

