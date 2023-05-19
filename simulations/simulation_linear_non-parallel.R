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

# c(1,2),c(1:3),c(1:4),c(1:5),c(1:10),
# c(1:20),c(1:35),c(1:55),c(1:100),c(1,6:15),
# c(1,6:50),c(1,2,6:15),c(1,2,6:50),c(1:3,6:15),c(1:3,6:50),
# c(1:4,6:15),c(1:4,6:50),c(6:20),c(6:30),c(6:45),
# c(6:70),c(6:100))

# loop settings
nsim <- 1000
n <- 100 # 1000
k <- c(2, seq(10, 100, 10)) # c(2, 10, 100, 500, 1000)
n_mod <- length(variables)

# create objects
results <- array(NA, dim = c(nsim, length(k), n_mod), dimnames = list(rep("sim", nsim), rep("k", length(k)), rep("mod", n_mod)))
mod <- list()
y_hat <- matrix(nrow = n, ncol = n_mod)

# loop 1
set.seed(2022)
for(j in 1:nsim){
  ind <- sample(1:N, n, replace = FALSE)
  for(l in 1:length(k)){ # choose lth fold number to test
    fold <- sample(1:k[l], n, replace = TRUE)
    for(i in 1:k[l]){ # outer cv loop
      for(m in 1:n_mod){print(c(j, l, i, m))
        mod[[m]] <- lm(y ~ ., cbind(x[ind, variables[[m]]][fold != i,], y = y[ind][fold != i]))
        y_hat[fold == i, m] <- predict(mod[[m]], x[ind, variables[[m]]][fold == i, ])
      }
    }
    for(m in 1:n_mod){
      results[j, l, m] <- mean((y[ind] - y_hat[, m])^2)
    }
  }
}

save(results, file = "/results20230223.RData")

# read results
theta.hat <- mse.hat <- readRDS("~/optimal_k_git/simulation_rnorm_lasso_new_results.rds")

# biassq + var = mse
biassq <- (apply(theta.hat, c(2, 3), mean) - theta)^2
var <- apply(theta.hat, c(2, 3), var)
mse <- apply((theta.hat - theta)^2, c(2, 3), mean)

model <- 3
plot(k, biassq[, model])
plot(k, var[, model])
plot(k, mse[, model])

results <- matrix(nrow = 5, ncol = 5)
for(sim in 1:5){print(sim)
  results[sim, ] <- apply(theta.hat[sim, , ], 1, which.min)
}
corrects <- colSums(results == 2)
plot(k, corrects)



