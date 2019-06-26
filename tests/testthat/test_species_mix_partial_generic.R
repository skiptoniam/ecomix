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

offset <- log(F6Data$Swept_Area)
spp_weights <- rep(1,ncol(y))
site_spp_weights <- matrix(1,nrow(y),ncol(y))
y_is_na <- matrix(FALSE,nrow(y),ncol(y))
G <- length(simulated_data$pi)
S <- length(simulated_data$sp.int)
nP <- ncol(X[,-1])
control <- species_mix.control()


library(ecomix)
load("~/Dropbox/ecomix_dev/developing_functions/kampala_code/Aug16.RData")
head(F6Data)
archetype_form <- as.formula(paste0('cbind(',paste(colnames(F6Data)[grep("spp",colnames(F6Data))][1:10],collapse = ','),")~cumTW.bs1+cumTW.bs2+cumTW.bs3"))
species_form <- ~ 1 + bathy1 + bathy2 + str1 + bathy_str + str2

test_dat <- make_mixture_data(F6Data[,grep("spp",colnames(F6Data))],F6Data[,(ncol(F6Data)-7):ncol(F6Data)])

# y <- as.matrix(F6Data[,grep("spp",colnames(F6Data))][1:50])
# X <- as.matrix(cbind(1,F6Data[,(ncol(F6Data)-7):(ncol(F6Data)-5)]))
# W <- as.matrix(F6Data[,(ncol(F6Data)-4):(ncol(F6Data))])

test_dat <- ecomix:::clean_data_sam(test_dat,archetype_form,species_form,distribution = 'negative_binomial')
test_dat$mf.X
test_dat$mf.W

y <- as.matrix(F6Data[,grep("spp",colnames(F6Data))])
spp_to_keep <- which(colSums(y)>10)
y <- y[,spp_to_keep]
X <- ecomix:::get_X_sam(archetype_form, test_dat$mf.X)
W <- ecomix:::get_W_sam(species_form, test_dat$mf.W)

n_mixtures <- G <- 4
S <- ncol(y)
spp_weights <- rep(1,S)
site_spp_weights <- matrix(1,nrow(y),S)
disty <- 4
y_is_na <- is.na(y)
inits <- NULL
control <- species_mix.control(quiet = FALSE)
offset <- log(F6Data$Swept_Area)




fm_sp_mods <-  surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W,
                                     site_spp_weights, offset, y_is_na, disty,
                                     .parallel = control$cores, .verbose = FALSE)


starting_values <- get_initial_values_sam(y = y, X = X, W = W,
                                                  spp_weights = spp_weights,
                                                  site_spp_weights = site_spp_weights,
                                                  offset = offset, y_is_na = y_is_na,
                                                  G = G, S = S,
                                                  disty = disty,
                                                  control = control)

fits <- starting_values$fits
taus <- starting_values$taus
pis <- starting_values$pis
first_fit <- starting_values$first_fit
logls_mus <- get_logls_partial_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)


ss <- 1
test <- apply_glm_spp_coefs_partial_sams(ss, y, X, W, G, taus,
                                         site_spp_weights,
                                         offset, y_is_na, disty, fits)
gg <- 1
test <- apply_glm_mix_coefs_partial_sams(gg, y, X, W, site_spp_weights,
                                         offset, y_is_na, disty, taus, fits, logls_mus$fitted)


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

