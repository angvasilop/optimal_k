
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
variables<-list(1:3,
                1:5,
                1:20,
                1:100,
                6:100)
                # c(1,2),c(1:3),c(1:4),c(1:5),c(1:10),
                # c(1:20),c(1:35),c(1:55),c(1:100),c(1,6:15),
                # c(1,6:50),c(1,2,6:15),c(1,2,6:50),c(1:3,6:15),c(1:3,6:50),
                # c(1:4,6:15),c(1:4,6:50),c(6:20),c(6:30),c(6:45),
                # c(6:70),c(6:100))

# loop settings
n <- 1000
nsim <- 1000
k <- c(2, 10, 100, 500, 1000)
n_mod <- length(variables)
results <- array(NA, dim = c(nsim, length(k), n_mod))
mod <- list()
y_hat <- matrix(nrow = n, ncol = n_mod)

# load glmnet
library(glmnet)

# create objects
best.lam <- c()
mod <- list()
predictions <- c()

# loop 1
set.seed(2022)
for(j in 1:nsim){
  ind <- sample(1:N, n, replace = FALSE)
  #dat <- as.matrix(data.frame(Y = y[ind],x[ind,]))
  for(l in 1:length(k)){ # choose lth fold number to test
    fold <- sample(1:k[l], n, replace = TRUE)
    for(i in 1:k[l]){ # outer cv loop
      for(m in 1:n_mod){print(c(j, l, i, m))
        best.lam[m] <- cv.glmnet(as.matrix(x[ind, variables[[m]]][fold != i,]), y[ind][fold != i], alpha = 1, folds = 5)$lambda.min # inner cv loop
        mod[[m]] <- glmnet(as.matrix(x[ind, variables[[m]]][fold != i,]), y[ind][fold != i], alpha = 1, lambda = best.lam[m])
        y_hat[fold == i,m] <- predict(mod[[m]], as.matrix(x[ind, variables[[m]]][fold == i,]))
      }
    }
    for(m in 1:n_mod){
      results[j, l, m] <- mean((y[ind] - y_hat[, m])^2)
    }
  }
}


save(results, file = "/results20230223.RData")











