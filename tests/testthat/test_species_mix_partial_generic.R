library(ecomix)
set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5), w2=rnorm(100,-1,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- species_mix.simulate(sam_form, ~w2, dat, theta, dist="bernoulli")

y <- simulated_data$species_data
X <- simulated_data$covariate_data
mf <- as.data.frame(cbind(y,X))
archetype_formula <- sam_form
species_formula <- ~ w2
distribution <- 'bernoulli'

offset <- rep(0,nrow(y))
spp_weights <- rep(1,ncol(y))
site_spp_weights <- matrix(1,nrow(y),ncol(y))
y_is_na <- matrix(FALSE,nrow(y),ncol(y))
G <- length(simulated_data$pi)
S <- length(simulated_data$sp.int)
nP <- ncol(X[,-1])
control <- species_mix.control()

head(F6Data)

archetype_form <- as.formula(paste0('cbind(',paste(colnames(F6Data)[grep("spp",colnames(F6Data))][1:10],collapse = ','),")~cumTW.bs1+cumTW.bs2+cumTW.bs3"))
species_form <- ~ 1 + bathy1 + bathy2 + str1 + bathy_str + str2

test_dat <- make_mixture_data(F6Data[,grep("spp",colnames(F6Data))],F6Data[,(ncol(F6Data)-7):ncol(F6Data)])

y <- as.matrix(F6Data[,grep("spp",colnames(F6Data))][1:100])
X <- as.matrix(F6Data[,(ncol(F6Data)-7):(ncol(F6Data)-5)])
W <- as.matrix(F6Data[,(ncol(F6Data)-4):(ncol(F6Data))])
n_mixtures <- G <- 5
S <- 100
spp_weights <- rep(1,S)
site_spp_weights <- matrix(1,nrow(y),S)
disty <- 4
y_is_na <- is.na(y)
inits <- NULL
control <- species_mix.control()

test <- fitmix_ECM_partial_sam(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)

tmp <- species_mix_partial.fit(y=y, X=X, W=W, G=n_mixtures, S=S,
                               spp_weights=spp_weights,
                               site_spp_weights=site_spp_weights,
                               offset=offset, disty=disty, y_is_na=y_is_na,
                               control=control, inits=inits)


test <- species_mix_partial(archetype_form, species_form, test_dat,
                            n_mixtures = 2, distribution="negative_binomial",
                            offset = NULL, weights = NULL,
                            bb_weights = NULL, control = NULL,
                            inits=NULL, standardise = FALSE,
                            titbits = TRUE)

