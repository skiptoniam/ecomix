require( RandomFields )
library(ggplot2)
location <- NULL
nsites <- 100
nspecies <- 50
n_mixtures <- 4
if(is.null(location)) locations <- cbind( "x"=runif(nsites, min=0,max=1), "y"=runif(nsites, min=0,max=1) )
dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
dat[,-1] <- scale(dat[,-1])
dat <- cbind(locations,dat)

## setup a matern spatial covariance for a spatial parameter
sigma = .01
kappa = 1
rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma^2, scale = 1 / kappa)
plot(rf_omega)
rf_sim <- function(model, x, y) {
  set.seed(sample.int(1e5L, 1L))
  suppressMessages(
    RandomFields::RFsimulate(model = model, x = x, y = y)$variable1
  )
}
omega_s <- rf_sim(model = rf_omega, locations[,'x'], locations[, "y"])
omega_s <- omega_s - mean(omega_s)
ggplot(data = data.frame(locations,omega_s),)+
  geom_point(aes(x,y,colour=omega_s))

# come up with a way to add a spatial formula.

archetype_formula <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
species_formula <- as.formula(~1+w1+w2)

beta <- matrix(c(-3.6,0.5,
                 -0.9,1.0,
                 0.9,-2.9,
                 2.2,5.4),
               4,2,byrow=TRUE)
gamma <- matrix(c(rnorm(50,1),rnorm(50)),50,2)

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

# # Putting spatial random effects on the species intercepts
# A_sp = matrix(NA, nrow=nsites, ncol=nspecies)
# sigma_sp <- 0.1
# kappa_sp <- 1
# model_A <- RandomFields::RMmatern(nu = 1, var = sigma_sp^2, scale = 1 / kappa_sp)
# for(ss in 1:S){
#   A_sp[,ss] = RandomFields::RFsimulate(model=model_A, x=locations[,'x'], y=locations[,'y'])@data[,1]
#   A_sp[,ss] = A_sp[,ss] - mean(A_sp[,ss]) + alpha[ss]
# }
# ggplot(data = data.frame(locations,A_sp),)+
#   geom_point(aes(x,y,colour=X4))

## simulate the groups and fitted values.
fitted <- matrix(0, dim(X)[1], S)
eta <- matrix(0, dim(X)[1], S)
group <- rep(0, S)
for (ss in seq_len(S)) {
  gg <- ceiling(stats::runif(1) * G)
  eta_spp <- W %*% c(alpha[ss],gamma[ss,])
  eta_mix <- X %*% beta[gg, ]
  eta[,ss] <- eta_spp + eta_mix + omega_s + offset
  fitted[, ss] <- link$linkinv(eta[,ss])
  group[ss] <- gg
}

matplot(eta)
points(omega_s,col="red",pch=16)

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

# adding in knots
# make_spde_w_knots <- function(x, y, n_knots, seed = 42, mesh = NULL) {
#   loc_xy <- cbind(x, y)
#
#   if (is.null(mesh)) {
#     if (n_knots >= nrow(loc_xy)) {
#       warning(
#         "Reducing `n_knots` to be one less than the ",
#         "number of data points."
#       )
#       n_knots <- nrow(loc_xy) - 1
#     }
#     set.seed(seed)
#     knots <- stats::kmeans(x = loc_xy, centers = n_knots)
#     loc_centers <- knots$centers
#     mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
#   } else {
#     knots <- list()
#     knots$cluster <- vapply(seq_len(nrow(loc_xy)), function(i)
#       RANN::nn2(mesh$loc[, 1:2, drop = FALSE],
#                 t(as.numeric(loc_xy[i, , drop = FALSE])),
#                 k = 1L
#       )$nn.idx,
#       FUN.VALUE = 1L
#     )
#     loc_centers <- NA
#   }
#   spde <- INLA::inla.spde2.matern(mesh)
#   A <- INLA::inla.spde.make.A(mesh, loc = loc_xy)
#   list(
#     x = x, y = y, mesh = mesh, spde = spde, cluster = knots$cluster,
#     loc_centers = loc_centers, A = A
#   )
# }
#
# plot_spde <- function(object) {
#   plot(object$mesh, main = NA, edge.color = "grey60", asp = 1)
#   points(object$x, object$y, pch = 21, col = "#00000070")
#   points(object$loc_centers, pch = 20, col = "red")
# }
#
# loc_xy <- cbind(locations[,"x"], locations[,"y"])
# bnd <- INLA::inla.nonconvex.hull(as.matrix(loc_xy), convex = -0.05)
# mesh <- INLA::inla.mesh.2d(
#   boundary = bnd,
#   max.edge = c(1, 1),
#   offset = -0.05,
#   # cutoff = c(2, 5),
#   # min.angle = 10
# )
# sp2 <- make_spde_w_knots(locations[,1],locations[,2], mesh = mesh)
# plot_spde(sp2)
# spde <- make_spde_w_knots(locations[,1],locations[,2],25)
# plot_spde(spde)

loc = locations
boundary = INLA::inla.nonconvex.hull(loc)
boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.35)
mesh = INLA::inla.mesh.2d(
  loc=loc,
  boundary = list(boundary,boundary2),
  max.edge=c(0.05, 0.2),
  cutoff=0.05
)
plot(mesh)
points(locations,col="tomato",pch=16)
A = INLA::inla.spde.make.A(mesh,loc)
spde = INLA::inla.spde2.matern(mesh, alpha=2)
spdeMatrices = spde$param.inla[c("M0","M1","M2")]

dats = list(Y = as.matrix(outcomes), # Response
            y_is_na = matrix(as.integer(!outcomes),nrow(X),S),
            X = X, # Design matrix for archetypes
            W = W, # Design matrix for species
            size = size,
            offy = offset, #offy is the offset indexed by sites (i)
            wts = site_spp_weights, #wts is a matrix indexed by sites, species (i,j).
            bb_wts = rep(1,S),
            nObs= nrow(X),# nsites.
            nG = G,       # n groups
            nS = S,
            family = as.integer(disty),
            link= as.integer(1),
            meshidxloc = mesh$idx$loc - 1,
            A = A,
            spdeMatrices = spdeMatrices)

gamma <- t(cbind(attr(simulated_data,"alpha"),attr(simulated_data,"gamma")))
beta[] <- beta[]+rnorm(length(beta),0,.2)
gamma[] <- gamma[]+rnorm(length(gamma),0,0.2)


#Define parameter object given to TMB-------
pars = list(beta=t(beta),
            gamma=gamma,
            eta = ecomix:::additive_logistic(attr(simulated_data,"pi"),TRUE)[-G],
            theta = rep(1,S),#1/exp(attr(simulated_data,'theta')))
            log_tau = log(0.01),
            log_kappa = log(1),
            x = rep(0.0, nrow(dats$spdeMatrices$M0)))


library(TMB)
compile("/home/woo457/Dropbox/ecomix/devsrc/tmbsam_gmrf.cpp","&> /tmp/logfile.log")
dyn.load(dynlib("/home/woo457/Dropbox/ecomix/devsrc/tmbsam_gmrf"))
obj <- MakeADFun(dats, pars, random="x", DLL="tmbsam_gmrf")
opt <- nlminb(obj$par, obj$fn, obj$gr)
SD = sdreport( obj )

g <- as.numeric(obj$gr(opt$par))
h <- stats::optimHess(opt$par, fn = obj$fn, gr = obj$gr)
tmb_opt$par <- tmb_opt$par - solve(h, g)
tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)

# Step 4 -- Simulate from predictive distribution
match_index = grep( "b_j", names(Opt$par) )
bhat_rj = mvtnorm::rmvnorm( n=1e4, mean=Opt$par[match_index], sigma=SD$cov.fixed[match_index,match_index] )

# predict response for new values
Xpred_z = seq( from=-10, to=10, length=1000 )
Ybounds_zj = matrix( NA, ncol=2, nrow=length(Xpred_z) )
for( z in 1:nrow(Ybounds_zj) ){
  ysim_r = bhat_rj[,1] + bhat_rj[,2]*Xpred_z[z]
  Ybounds_zj[z,] = quantile( ysim_r, prob=c(0.1,0.9) )
}

# plot results
plot( x=X, y=Y )
abline( a=Opt$par[match_index][1], b=Opt$par[match_index][2] )
polygon( x=c(Xpred_z,rev(Xpred_z)), y=c(Ybounds_zj[,1],rev(Ybounds_zj[,2])), col=rgb(1,0,0,0.2) )


