context('species_mix poisson')

library(ecomix)


testthat::test_that('species mix poisson', {

set.seed(42)
form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
theta <- matrix(c(-2.9,-1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- simulate_species_mix_data(form,~1,dat,theta,dist="poisson")

y <- simulated_data$species_data
X <- simulated_data$covariate_data
# X[,-1] <- scale(X[,-1])
offset <- rep(0,nrow(y))
weights <- rep(1,nrow(y))
spp_wts <- rep(1,ncol(y))
site_spp_wts <- matrix(1,nrow(y),ncol(y))
y_is_na <- matrix(FALSE,nrow(y),ncol(y))
G <- length(simulated_data$pi)
S <- length(simulated_data$sp.int)

i <- 1
testthat::expect_length(ecomix:::apply_glmnet_poisson(i, y, X, weights, offset),3)


## get to this in a bit
gg <- 1
ecomix:::apply_glm_poisson_group_tau(gg, y, X, y_is_na, tau)


})
