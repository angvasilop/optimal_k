
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
n <- 500 # 1000
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
set.seed(2022)
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
out <- readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_500_n_500.rds")
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
df <- data.frame(k, corrects)
plot(df$k, df$corrects)

# nls
nls.fit <- nls(corrects ~ a*k^-1 + b, start = list(a = 14.5, b = 660))
.a <- coef(nls.fit)[1]
.b <- coef(nls.fit)[2]
eq <- function(x){.a*(x)^-1 + .b}
lines(eq(0:700), type = "l")

# Mann-Kendall tests
df$pred <- predict(nls.fit, df$k)
df$pval <- c()
library(Kendall)
for(i in 1:(length(df$pred) - 2)){
  df$pval[i] <- MannKendall(df$pred[i:length(df$pred)])$sl
}
df$pval[(length(df$pred) - 1):length(df$k)] <- 0
length(df$pval)
df$k[which.max(df$pval[df$pval < 0.05])]
df$k[which.max(df$pval[df$pval < 0.05])]/n

# k vs. n
sample.size <- seq(100, 700, 100)
k.opt.corrected <- c(0.4,0.65,0.767,0.825,0.86, 0.867,0.886)
plot(sample.size, k.opt.corrected, xlab = "sample.size (n)", ylab = "k.optimal/n")

# second fit
nls.fit <- nls(k.opt.corrected ~ a*sample.size^-1 + b, start = list(a = 1, b = 0))
.a <- coef(nls.fit)[1]
.b <- coef(nls.fit)[2]
eq <- function(x){.a*(x)^-1 + .b}
lines(eq(0:n), type = "l")

