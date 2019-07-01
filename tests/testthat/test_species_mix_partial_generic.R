library(ecomix)
set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~x1+x2"))
spp_form <- as.formula(~1+w2)
beta <- matrix(c(-3.6,0.5,
                 -0.9,1.0,
                  0.9,-2.9,
                  2.2,5.4),
                4,2,byrow=TRUE)
gamma <- rnorm(50,1)
dat <- data.frame(y=rep(1,100), x1=runif(100,0,2.5), x2=rnorm(100,0,2.5), w2=rnorm(100,-1,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- species_mix.simulate(sam_form, spp_form, dat = dat,
                                       beta = beta, gamma = gamma,
                                       n_mixtures = 4,
                                       distribution = "bernoulli")

archetype_formula <- sam_form
species_formula <- spp_form

test_dat <- ecomix:::clean_data_sam(simulated_data, archetype_formula, species_formula, distribution = 'bernoulli')
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
y_is_na <- is.na(y)
inits <- NULL
control <- species_mix.control(quiet = FALSE)
offset <- rep(0,nrow(X))

fm_sp_mods <-  surveillance::plapply(seq_len(S), ecomix:::apply_glmnet_sam_inits, y, X, W,
                                     site_spp_weights, offset, y_is_na, disty,
                                     .parallel = control$cores, .verbose = FALSE)


starting_values <- ecomix:::get_initial_values_sam(y = y, X = X, W = W,
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
logls_mus <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)


ss <- 1
spp_conditional_max <- ecomix:::apply_glm_spp_coefs_sams(ss, y, X, W, G, taus,
                                                         site_spp_weights,
                                                         offset, y_is_na, disty, fits)
gg <- 1
mix_conditional_max <- ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W, site_spp_weights,
                                             offset, y_is_na, disty, taus, fits, logls_mus$fitted)

partial_ECM <- ecomix:::fitmix_ECM_sam(y, X, W, spp_weights, site_spp_weights,
                                       offset, y_is_na, G, S, disty,
                                       control=species_mix.control(em_steps=10))

start_vals <- ecomix:::get_starting_values_sam(y = y, X = X, W = W,
                                      spp_weights = spp_weights,
                                      site_spp_weights = site_spp_weights,
                                      offset = offset,
                                      y_is_na = y_is_na,
                                      G = G, S = S,
                                      disty = disty,
                                      control = control)


tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,
                             S,G,disty,start_vals,
                             control=ecomix:::species_mix.control(optimise_cpp = 0,
                                                         loglOnly_cpp = 1))




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
spp_to_keep <- which(colSums(y>1)>20)
y <- y[,spp_to_keep]
X <- ecomix:::get_X_sam(archetype_form, test_dat$mf.X)
W <- ecomix:::get_W_sam(species_form, test_dat$mf.W)

n_mixtures <- G <- 6
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

starting_values <- ecomix:::get_initial_values_sam(y = y, X = X, W = W,
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
logls_mus <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty, get_fitted = TRUE)


ss <- 1
test <- ecomix:::apply_glm_spp_coefs_sams(ss, y, X, W, G, taus,
                                         site_spp_weights,
                                         offset, y_is_na, disty, fits)
gg <- 1
test <- ecomix:::apply_glm_mix_coefs_sams(gg, y, X, W, site_spp_weights,
                                         offset, y_is_na, disty, taus, fits, logls_mus$fitted)

test <- fitmix_ECM_sam(y, X, W, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)

start_vals <- ecomix:::get_starting_values_sam(y = y, X = X, W = W,
                                               spp_weights = spp_weights,
                                               site_spp_weights = site_spp_weights,
                                               offset = offset,
                                               y_is_na = y_is_na,
                                               G = G, S = S,
                                               disty = disty,
                                               control = control)

tmp <- ecomix:::sam_optimise(y,X,W,offset,spp_weights,site_spp_weights,y_is_na,S,G,disty,start_vals,control)



tmp <- species_mix.fit(y=y, X=X, W=W, G=n_mixtures, S=S,
                               spp_weights=spp_weights,
                               site_spp_weights=site_spp_weights,
                               offset=offset, disty=disty, y_is_na=y_is_na,
                               control=control, inits=inits)


test <- species_mix(archetype_form, species_form, test_dat,
                            n_mixtures = 2, distribution="negative_binomial",
                            offset = NULL, weights = NULL,
                            bb_weights = NULL, control = NULL,
                            inits=NULL, standardise = FALSE,
                            titbits = TRUE)

