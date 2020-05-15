library('mgcv')
set.seed(20)
n <- 100
S <- 20
beta <- 4

alpha <- runif(S,-7,-4) #intercept
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
f3 <- function(x) 0 * x
theta <-  log( 1 + rgamma( n=S, shape=1, scale=0.75))

x0 <- runif(n)
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)

X <- cbind(x0,x1,x2,x3)

fb1 <- f0(x0) + f1(x1) + f2(x2)
fb2 <- f0(x0) + f2(x1)
fb3 <- f2(x0) + f1(x1)
fb4 <- f1(x1) + f3(x3)
scale <- 1

f <- cbind(fb1,fb2,fb3,fb4)
fitted <- matrix(0, dim(X)[1], S)
group <- rep(0, S)
offset <- rep(0,dim(X)[1])
link <- make.link("log")
for (ss in seq_len(S)) {
  gg <- ceiling(stats::runif(1) * G)
  eta_spp <- alpha[ss]
  eta_mix <- f[,gg]
  eta <- eta_spp + eta_mix + offset
  fitted[, ss] <- link$linkinv(eta*scale)
  group[ss] <- gg
}

outcomes <- matrix(rnbinom(n * S, mu=as.numeric( fitted), size=1/rep(exp(theta), each=n)), nrow = n, ncol = S)
pi <- tapply(group, group, length)/S

Y <- outcomes
X <- X
W <- matrix(1,nrow(X),1)
offy <- offset
spp_site_wts <- matrix(1,nrow(X),S)
spp_wts <- rep(1,S)

df <- data.frame(y = Y[,2],X)

group

m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = df)
set.seed(10)
nsim <- 50
drawx <- function(x, n) runif(n, min = min(x), max = max(x))
newx <- data.frame(apply(X,2,drawx, n = nsim))

sig2 <- m$sig2
set.seed(10)

newx <- transform(newx, newy = predict(m, newx, type = "response"))
newx <- transform(newx, ysim = rnorm(nsim, mean = newy, sd = sqrt(sig2)))

pred <- data.frame(apply(df[,-1],2, function(x) seq(min(x), max(x), length = 500)))
pred <- transform(pred, fitted = predict(m, newdata = pred,
                                         type = "response"))
library('ggplot2')
theme_set(theme_bw())

ggplot(df) +
  geom_point(aes(x = x1, y = y), colour = "grey") +
  # geom_line(aes(x = x2, y = f), colour = "forestgreen", size = 1.3) +
  geom_line(aes(x = x1, y = fitted), data = pred, colour = "blue") +
  geom_point(aes(x = x1, y = newy), data = newx, colour = "blue", size = 2) +
  geom_point(aes(x = x1, y = ysim), data = newx, colour = "red",
             size = 2)
