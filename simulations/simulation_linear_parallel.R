
# create population
N <- 500000
x<-data.frame(matrix(nrow=500000,ncol=100))
set.seed(2022)
n_var <- 1:100
for(i in n_var){
  x[,i]<-rnorm(n=N,mean=0,sd=1)
}

# define parameters
beta0 <- 0
beta1 <- 1
beta2 <- 1
beta3 <- 1
beta4 <- 1
beta5 <- 1
sigma <- 10

# define theta
Ey <- beta0 + beta1*x[,1] + beta2*x[,2] + beta3*x[,3] + beta4*x[,4] + beta5*x[,5]
y <- Ey + rnorm(nrow(x),0,sigma)
theta <- mean((y - Ey)^2)

# choose models
variables<-list(1:3, 1:5, 1:20, 1:100, 6:100)

# loop settings
n <- 100
nsim <- 1000
k <- c(2, seq(10, n, 10)) # c(2, 10, 100, 500, 1000)
n_mod <- length(variables)

# create objects
best.lam <- c()
mod <- list()
y_hat <- matrix(nrow = n, ncol = n_mod)
results <- matrix(nrow = length(k), ncol = n_mod)
out <- list()

# loop
library(foreach)
library(doParallel)

# register cluster
registerDoParallel(makeCluster(detectCores() - 1))

# loop
set.seed(2025)
out <- foreach(j = 1:nsim, .combine = "list") %dopar% {
  ind <- sample(1:N, n, replace = FALSE)
  for(l in 1:length(k)){
    fold <- sample(1:k[l], n, replace = TRUE)
    for(i in 1:k[l]){print(c(j, l, i))
      for(m in 1:n_mod){
        mod <- lm(y ~ ., cbind(x[ind, variables[[m]]][fold != i,], y = y[ind][fold != i]))
        y_hat[fold == i, m] <- predict(mod, x[ind, variables[[m]]][fold == i, ])
      }
    }
    for(m in 1:n_mod){
      results[l, m] <- mean((y[ind] - y_hat[, m])^2)
    }
  }
  return(results)
}

# read results
out <- readRDS("~/Downloads/parallel_simulation_linear_output_k_2_100_n_100_d.rds")
theta.hat <- mse.hat <- array(unlist(out), dim = c(length(k), n_mod, nsim), dimnames = list(rep("k", length(k)), rep("mod", n_mod), rep("sim", nsim)))

# biassq + var = mse
biassq <- (apply(theta.hat, c(1, 2), mean) - theta)^2
var <- apply(theta.hat, c(1, 2), var)
mse <- apply((theta.hat - theta)^2, c(1, 2), mean)

# plots
library(plotrix)

# bias^2 plots
plot(k, biassq[, 1], type = "b", ylim = c(min(biassq), max(biassq)), ylab = "bias^2")
points(k, biassq[, 2], type = "b", col = "blue")
points(k, biassq[, 3], type = "b", col = "red")
points(k, biassq[, 4], type = "b", col = "purple")
points(k, biassq[, 5], type = "b", col = "green")
legend(600, 250, legend=c("1:3", "1:5", "1:20", "1:100", "6:100"),
       col=c("black", "blue", "red", "purple", "green"), lty = 1, cex = 0.8)

# variance plot
plot(k, var[, 1], type = "b", ylim = c(min(var), max(var)), ylab = "var")
points(k, var[, 2], type = "b", col = "blue")
points(k, var[, 3], type = "b", col = "red")
points(k, var[, 4], type = "b", col = "purple")
points(k, var[, 5], type = "b", col = "green")
legend(600, 50, legend=c("1:3", "1:5", "1:20", "1:100", "6:100"),
       col=c("black", "blue", "red", "purple", "green"), lty = 1, cex = 0.8)

# MSE plot
plot(k, mse[, 1], type = "b", ylim = c(min(mse), max(mse)), ylab = "MSE")
points(k, mse[, 2], type = "b", col = "blue")
points(k, mse[, 3], type = "b", col = "red")
points(k, mse[, 4], type = "b", col = "purple")
points(k, mse[, 5], type = "b", col = "green")
legend(600, 250, legend=c("1:3", "1:5", "1:20", "1:100", "6:100"),
       col=c("black", "blue", "red", "purple", "green"), lty = 1, cex = 0.8)

# which k results in the correct model having the lowest MSE the most often?
results <- matrix(nrow = nsim, ncol = length(k))
for(sim in 1:nsim){
  results[sim, ] <- apply(theta.hat[,, sim], 1, which.min)
}
corrects <- colSums(results == 2)
plot(k, corrects)

# nls
nls.fit <- nls(corrects ~ c + d*k^-1, start = list(c = 970, d = -208))
.c <- coef(nls.fit)[1]
.d <- coef(nls.fit)[2]
eq <- function(x){.c + .d*x^-1}
pred <- eq(2:n)
lines(pred, type = "l")

# elbow point, optimal.k/n vs. n
library("smerc")
elbow_point(2:n, pred)$x
samp.n <- seq(100, 1000, 100)
opt.k <- c(14, 20, 24, 28, 32, 35, 37, 40, 42, 45) # departing from LOOCV
plot(samp.n, opt.k/samp.n)



opt.k.elb <- c(20/100, 20/200, 30/300, 30/400, 30/500, 30/600, 30/700, 20/800, 10/900)
opt.k.lo <- c(79/100, 52/200, 75/300, 104/400, 128/500, 159/600, 180/700, 208/800, 233/900)

# Mann-Kendall tests
df <- data.frame(k, corrects, pred)
df$pval <- c()
library(Kendall)
for(i in 1:(length(df$pred) - 2)){
  df$pval[i] <- MannKendall(df$pred[i:length(df$pred)])$sl
}
df$pval[(length(df$pred) - 1):length(df$k)] <- 0
df$k[which.max(df$pval[df$pval < 0.05/(n/10 - 1)])]
df$k[which.max(df$pval[df$pval < 0.05/(n/10 - 1)])]/n
# k vs. n
sample.size <- seq(100, 900, 100)
k.opt.corrected <- c(0.4,0.65,0.767,0.825,0.86, 0.883, 0.886, 0.9, 0.911)
plot(sample.size, k.opt.corrected, xlab = "sample.size (n)", ylab = "k.optimal/n")
logit(k.optimal.corrected)

ss.fit <- smooth.spline(k, corrects)
pred <- predict(ss.fit, 2:n)$y
lines(pred)
elbow_point(2:n, pred)$x
samp.n <- seq(200, 900, 100)
opt.k.ss <- c(14, 24, 48, 28, 31, 6, 16, 26)
