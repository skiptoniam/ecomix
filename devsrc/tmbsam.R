## test tmb sam.
library(ecomix)
set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
spp_form <- as.formula(~1+w1+w2)
beta <- matrix(c(-3.6,0.5,
                 -0.9,1.0,
                 0.9,-2.9,
                 2.2,5.4),
               4,2,byrow=TRUE)
gamma <- matrix(c(rnorm(50,1),rnorm(50,-2)),50,2)
dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5),w1=rnorm(100,2,1), w2=rnorm(100,-1,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- species_mix.simulate(sam_form, spp_form, dat = dat,
                                       beta = beta, gamma = gamma,
                                       n_mixtures = 4,
                                       distribution = "bernoulli")
attr(simulated_data,"pi")



archetype_formula <- sam_form
species_formula <- spp_form

test_dat <- ecomix:::clean_data_sam(simulated_data, archetype_formula, species_formula, distribution = 'negative_binomial')
test_dat$mf.X
test_dat$mf.W

y <- simulated_data[,1:50]
X <- ecomix:::get_X_sam(archetype_formula, test_dat$mf.X)
W <- ecomix:::get_W_sam(species_formula, test_dat$mf.W)

n_mixtures <- G <- 4
S <- ncol(y)
spp_weights <- rep(1,S)
site_spp_weights <- matrix(1,nrow(y),S)
disty <- 1
link <- 2
y_is_na <- is.na(y)
inits <- NULL
control <- species_mix.control(quiet = FALSE)
offset <- rep(0,nrow(X))
size <- as.integer(rep(1,nrow(X)))

library(TMB)
# compile("/home/woo457/Dropbox/ecomix/devsrc/tmbsam.cpp",flags="-O0 -g",DLLFLAGS="")
compile("/home/woo457/Dropbox/ecomix/devsrc/tmbsam.cpp","&> /tmp/logfile.log")
dyn.load(dynlib("/home/woo457/Dropbox/ecomix/devsrc/tmbsam"))

#Define data object which is given to TMB---
dats = list(Y = as.matrix(y), # Response
            y_is_na = matrix(as.integer(!y_is_na),nrow(X),S),
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
            link=as.integer(link),
            keep_mu = as.integer(0))#,

#-------------------------------------------

gamma <- t(cbind(attr(simulated_data,"alpha"),attr(simulated_data,"gamma")))
beta[] <- beta[]+rnorm(length(beta),0,.2)
gamma[] <- gamma[]+rnorm(length(gamma),0,0.2)


#Define parameter object given to TMB-------
pars = list(beta=t(beta),
           gamma=gamma,
           eta = ecomix:::additive_logistic(attr(simulated_data,"pi"),TRUE)[-G],
           theta = rep(1,S))#1/exp(attr(simulated_data,'theta')))

#   beta = rep(0,sum(Sdims)),  # Spline coefficients
#   log_lambda = rep(rep(0,length(Sdims))), #Log spline penalization coefficients
#   log_sigma = 0
# )
#-------------------------------------------

#Fit model----------------------------------
obj = MakeADFun(data = dats, parameters = pars, DLL = "tmbsam")
Opt = nlminb( start=obj$par, objective=obj$fn, gr=obj$gr)
SD = sdreport( obj )

# ### Predict the confedence intervals.
#
#
# # Generate data
# n_data = 30
# X = runif( n_data )
# Y = 1 + 0.2*X + rnorm( n_data )
#
# # Step 1 -- make and compile template file
# compile( "devsrc/TMB_intervals.cpp" )
#
# # Step 2 -- build inputs and object
# dyn.load( dynlib("devsrc/TMB_intervals") )
# Data = list( "y_i"=Y, "X_ij"=cbind(1,X) )
# Params = list( "b_j"=rep(0,ncol(Data$X_ij)), "log_sd"=0 )
# Obj = MakeADFun( data=Data, parameters=Params, DLL="TMB_intervals")
#
# # Step 3 -- test and optimize
# Opt = nlminb( start=Obj$par, objective=Obj$fn, gr=Obj$gr )
# SD = sdreport( Obj )
#
# # Step 4 -- Simulate from predictive distribution
# match_index = grep( "b_j", names(Opt$par) )
# bhat_rj = mvtnorm::rmvnorm( n=1e4, mean=Opt$par[match_index], sigma=SD$cov.fixed[match_index,match_index] )
#
# # predict response for new values
# Xpred_z = seq( from=-10, to=10, length=1000 )
# Ybounds_zj = matrix( NA, ncol=2, nrow=length(Xpred_z) )
# for( z in 1:nrow(Ybounds_zj) ){
#   ysim_r = bhat_rj[,1] + bhat_rj[,2]*Xpred_z[z]
#   Ybounds_zj[z,] = quantile( ysim_r, prob=c(0.1,0.9) )
# }
#
# # plot results
# plot( x=X, y=Y )
# abline( a=Opt$par[match_index][1], b=Opt$par[match_index][2] )
# polygon( x=c(Xpred_z,rev(Xpred_z)), y=c(Ybounds_zj[,1],rev(Ybounds_zj[,2])), col=rgb(1,0,0,0.2) )
#
