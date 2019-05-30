library(ecomix)
set.seed(42)
sam_form <- as.formula(paste0('cbind(',paste(paste0('spp',1:50),collapse = ','),")~1+x1+x2"))
theta <- matrix(c(-2.9,1.6,0.5,1,-0.9,1,.9,2.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x1=runif(100,0,2.5),x2=rnorm(100,0,2.5),w2=rnorm(100,-1,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- species_mix.simulate(sam_form,~1+w2,dat,theta,dist="bernoulli")

y <- simulated_data$species_data
X <- simulated_data$covariate_data
mf <- as.data.frame(cbind(y,X))
archetype_formula <- sam_form
species_formula <- ~ w2
distribution <- 'bernoulli'

offset <- rep(0,nrow(y))
# weights <- rep(1,nrow(y))
spp_weights <- rep(1,ncol(y))
site_spp_weights <- matrix(1,nrow(y),ncol(y))
y_is_na <- matrix(FALSE,nrow(y),ncol(y))
G <- length(simulated_data$pi)
S <- length(simulated_data$sp.int)
nP <- ncol(X[,-1])
control <- species_mix.control()

control <- species_mix_partial.control()
