
# load libraries
library(ggplot2)
library(ggpubr)
library(kableExtra)

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
n <- c(100, 250, 500, 1000)
nsim <- 1000
k <- list(n100 = seq(2, 100, 2),
          n250 = c(2, seq(10, 250, 10)),
          n500 = c(2, seq(10, 500, 10)),
          n1000 = c(2, 10, 100, 500, 1000)) # for n = 1000
n_mod <- length(variables)
out <- list(n100 = readRDS("~/optimal_k_git/results/lasso/parallel_simulation_lasso_output_k_2_100_n_100.rds"),
            n250 = readRDS("~/optimal_k_git/results/lasso/parallel_simulation_lasso_output_k_2_250_n_250.rds"),
            n500 = readRDS("~/optimal_k_git/results/lasso/parallel_simulation_lasso_output_k_2_500_n_500.rds"),
            n1000 = readRDS("~/optimal_k_git/results/lasso/parallel_simulation_lasso_output_k_2_1000_n_1000.rds"))
titles <- list("n = 100", "n = 250", "n = 500", "n = 1000")
plotsout <- list(list(), list(), list())
figures <- list()

for(i in 1:length(out)){
  # read results
  theta.hat <- mse.hat <- array(unlist(out[[i]]), dim = c(length(k[[i]]), n_mod, nsim))
  
  # biassq + var = mse
  biassq <- as.data.frame((apply(theta.hat, c(1, 2), mean) - theta)^2)
  var <- as.data.frame(apply(theta.hat, c(1, 2), var))
  mse <- as.data.frame(apply((theta.hat - theta)^2, c(1, 2), mean))
  metrics <- list(biassq, var, mse)
  
  for(j in 1:length(metrics)){
    # plot
    plotsout[[j]][[i]] <- ggplot(data.frame(k = k[[i]], metrics[[j]]), aes(x = k))+
      geom_line(aes(y = V1, color = "1 - 3")) +
      geom_line(aes(y = V2, color = "1 - 5")) +
      geom_line(aes(y = V3, color = "1 - 20")) +
      geom_line(aes(y = V4, color = "1 - 100")) +
      geom_line(aes(y = V5, color = "6 - 100")) +
      ggtitle(paste(titles[i])) +
      labs(x = NULL, y = NULL) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            axis.line=element_line(colour="black",size=0.4),
            axis.text=element_text(size=10,color="black"),
            plot.title=element_text(hjust=0.5,size=10, color = "black"),
            legend.position = c(1, 1),
            legend.justification = c("right", "top"),
            legend.key=element_rect(fill = NA),
            legend.title = element_text(size = 10, color = "black"),
            legend.text = element_text(size = 10, color = "black")) +
      scale_color_discrete(name = "Regression on variables", breaks=c("1 - 3", "1 - 5", "1 - 20", "1 - 100", "6 - 100"))
  }
}

for(i in 1:length(metrics)){
  figures[[i]] <- ggarrange(plotsout[[i]][[1]] + rremove("ylab") + rremove("xlab"),
                            plotsout[[i]][[2]] + rremove("ylab") + rremove("xlab"),
                            plotsout[[i]][[3]] + rremove("ylab") + rremove("xlab"),
                            plotsout[[i]][[4]] + rremove("ylab") + rremove("xlab"),
                            labels = NULL,
                            ncol = 2, nrow = 2,
                            common.legend = TRUE,
                            legend = "right",
                            align = "hv", 
                            font.label = list(size = 10, color = "black", family = NULL, position = "top"))
}



figure1 <- annotate_figure(figures[[1]],
                top = text_grob("Squared bias vs. k", size = 10, x = 0.4),
                bottom = text_grob("k", size = 10, x = 0.4),
                left = text_grob("Squared bias", size = 10, rot = 90))
figure2 <- annotate_figure(figures[[2]],
                top = text_grob("Variance vs. k", size = 10, x = 0.4),
                bottom = text_grob("k", size = 10, x = 0.4),
                left = text_grob("Variance", size = 10, rot = 90))
figure3 <- annotate_figure(figures[[3]],
                top = text_grob("Mean squared error (MSE) vs. k", size = 10, x = 0.4),
                bottom = text_grob("k", size = 10, x = 0.4),
                left = text_grob("MSE", size = 10, rot = 90))

counts <- list()
corrects <- c()
optplotsout <- list()
for(i in 1:length(n)){
  
  # read results
  theta.hat <- mse.hat <- array(unlist(out[[i]]), dim = c(length(k[[i]]), n_mod, nsim))
  results <- matrix(nrow = nsim, ncol = length(k[[i]]))
  for(sim in 1:nsim){
    results[sim, ] <- apply(theta.hat[,, sim], 1, which.min)
  }
  corrects <- counts[[i]] <- colSums(results == 2)
  
  # plot
  optplotsout[[i]] <- ggplot(data.frame(k = k[[i]], corrects), aes(x = k))+
    geom_line(aes(y = corrects)) +
    ggtitle(paste(titles[i])) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(colour="black",size=0.4),
          axis.text=element_text(size=10,color="black"),
          plot.title=element_text(hjust=0.5,size=10, color = "black"))
}

optfigures <- ggarrange(optplotsout[[1]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[2]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[3]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[4]] + rremove("ylab") + rremove("xlab"),
                     labels = NULL,
                     ncol = 2, nrow = 2,
                     common.legend = TRUE,
                     legend = "right",
                     align = "hv", 
                     font.label = list(size = 10, color = "black", family = NULL, position = "top"))

figure4 <- annotate_figure(optfigures,
                top = text_grob("Lowest error model selection", size = 10),
                bottom = text_grob("k", size = 10),
                left = text_grob("Number of times true model has lowest MSE", size = 10, rot = 90))

dffortable <- data.frame(n = c(100, 250, 500, 1000),
                         Minimum = c(min(counts[[1]]), min(counts[[2]]), min(counts[[3]]), min(counts[[4]])),
                         Maximum = c(max(counts[[1]]), max(counts[[2]]), max(counts[[3]]), max(counts[[4]])),
                         Range = c(max(counts[[1]]) - min(counts[[1]]), max(counts[[2]]) - min(counts[[2]]), max(counts[[3]]) - min(counts[[3]]), max(counts[[4]]) - min(counts[[4]])),
                         Variance = c(var(counts[[1]]), var(counts[[2]]), var(counts[[3]]), var(counts[[4]])))
