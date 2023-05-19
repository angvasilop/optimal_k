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
n <- 750 # 250 # 500 # 100 # 1000
nsim <- 1000
# k <- c(seq(2, 100, 2), 500, 1000)
# k <- seq(2, 100, 2) # for n = 100 and n = 1000
# k <- c(2, seq(10, 250, 10)) # for n = 250
k <- c(2, seq(10, 500, 10)) # for n = 500
# k <- c(2, seq(10, 750, 10)) # for n = 750
# k <- c(2, 10, 100, 500, 1000) # for n = 1000
n_mod <- length(variables)

# create objects
best.lam <- c()
mod <- list()
y_hat <- matrix(nrow = n, ncol = n_mod)
results <- matrix(nrow = length(k), ncol = n_mod)
out <- list()

# loop
library(foreach)
library(glmnet)
library(doParallel)

# register cluster
registerDoParallel(makeCluster(detectCores() - 1))

# loop
set.seed(2022)
out <- foreach(j = 1:nsim, .combine = "list", .packages = "glmnet") %dopar% {
  ind <- sample(1:N, n, replace = FALSE)
  for(l in 1:length(k)){
    fold <- sample(1:k[l], n, replace = TRUE)
    for(i in 1:k[l]){
      for(m in 1:n_mod){print(c(j, l, i, m))
        best.lam[m] <- cv.glmnet(as.matrix(x[ind, variables[[m]]][fold != i,]), y[ind][fold != i], alpha = 1, folds = 5)$lambda.min
        mod[[m]] <- glmnet(as.matrix(x[ind, variables[[m]]][fold != i,]), y[ind][fold != i], alpha = 1, lambda = best.lam[m])
        y_hat[fold == i, m] <- predict(mod[[m]], as.matrix(x[ind, variables[[m]]][fold == i, ]))
      }
    }
    for(m in 1:n_mod){
      results[l, m] <- mean((y[ind] - y_hat[, m])^2)
    }
  }
  return(results)
}

# read results
out <- readRDS("~/optimal_k_git/results/parallel_simulation_lasso_output_k_2_500_n_500 copy.rds")
theta.hat <- mse.hat <- array(unlist(out), dim = c(length(k), n_mod, nsim))

# biassq + var = mse
biassq <- (apply(theta.hat, c(1, 2), mean) - theta)^2
var <- apply(theta.hat, c(1, 2), var)
mse <- apply((theta.hat - theta)^2, c(1, 2), mean)

# plots
library(plotrix)

# bias^2 plot
plot(k, biassq[, 1], type = "b", ylim = c(min(biassq), max(biassq)), ylab = "bias^2")
points(k, biassq[, 2], type = "b", col = "blue")
points(k, biassq[, 3], type = "b", col = "red")
points(k, biassq[, 4], type = "b", col = "purple")
points(k, biassq[, 5], type = "b", col = "green")
legend(800, 20, legend=c("1:3", "1:5", "1:20", "1:100", "6:100"),
       col=c("black", "blue", "red", "purple", "green"), lty = 1, cex = 0.8)

# variance plot
plot(k, var[, 1], type = "b", ylim = c(min(var), max(var)), ylab = "var")
points(k, var[, 2], type = "b", col = "blue")
points(k, var[, 3], type = "b", col = "red")
points(k, var[, 4], type = "b", col = "purple")
points(k, var[, 5], type = "b", col = "green")
legend(800, 22, legend=c("1:3", "1:5", "1:20", "1:100", "6:100"),
       col=c("black", "blue", "red", "purple", "green"), lty = 1, cex = 0.8)

# MSE plot
plot(k, mse[, 1], type = "b", ylim = c(min(mse), max(mse)), ylab = "MSE")
points(k, mse[, 2], type = "b", col = "blue")
points(k, mse[, 3], type = "b", col = "red")
points(k, mse[, 4], type = "b", col = "purple")
points(k, mse[, 5], type = "b", col = "green")
legend(800, 30, legend=c("1:3", "1:5", "1:20", "1:100", "6:100"),
       col=c("black", "blue", "red", "purple", "green"), lty = 1, cex = 0.8)

# which k results in the correct model having the lowest MSE the most often?
results <- matrix(nrow = nsim, ncol = length(k))
for(sim in 1:nsim){
  results[sim, ] <- apply(theta.hat[,, sim], 1, which.min)
}
corrects <- colSums(results == 2)
plot(k, corrects, type = "b", ylab = "How many times correct model had lowest MSE (n = 250)")
df <- data.frame(k, corrects)

# plot
samplesize <- c(100, 250, 500, 1000)
optimalk <- c(24, 80, 70, 4)
abline(plot(samplesize, optimalk))


