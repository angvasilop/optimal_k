
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
variables <- list(1:3, 1:5, 1:20, 1:100, 6:100)

# loop settings
n <- seq(100, 1000, 100)
nsim <- 1000
k <- list(n100 = c(2, seq(10, 100, 10)),
          n200 = c(2, seq(10, 200, 10)),
          n300 = c(2, seq(10, 300, 10)),
          n400 = c(2, seq(10, 400, 10)),
          n500 = c(2, seq(10, 500, 10)),
          n600 = c(2, seq(10, 600, 10)),
          n700 = c(2, seq(10, 700, 10)),
          n800 = c(2, seq(10, 800, 10)),
          n900 = c(2, seq(10, 900, 10)),
          n1000 = c(2, seq(10, 1000, 10)))
n_mod <- length(variables)
out <- list(n100 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_100_n_100.rds"),
            n200 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_200_n_200.rds"),
            n300 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_300_n_300.rds"),
            n400 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_400_n_400.rds"),
            n500 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_500_n_500.rds"),
            n600 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_600_n_600.rds"),
            n700 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_700_n_700.rds"),
            n800 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_800_n_800.rds"),
            n900 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_900_n_900.rds"),
            n1000 = readRDS("~/optimal_k_git/results/linear/parallel_simulation_linear_output_k_2_1000_n_1000.rds"))
titles <- list("n = 100", "n = 200", "n = 300", "n = 400", "n = 500", "n = 600", "n = 700", "n = 800", "n = 900", "n = 1000")
plotsout <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
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
                            plotsout[[i]][[5]] + rremove("ylab") + rremove("xlab"),
                            plotsout[[i]][[6]] + rremove("ylab") + rremove("xlab"),
                            plotsout[[i]][[7]] + rremove("ylab") + rremove("xlab"),
                            plotsout[[i]][[8]] + rremove("ylab") + rremove("xlab"),
                            plotsout[[i]][[9]] + rremove("ylab") + rremove("xlab"),
                            plotsout[[i]][[10]] + rremove("ylab") + rremove("xlab"),
                            labels = NULL,
                            ncol = 2, nrow = 5,
                            common.legend = TRUE,
                            legend = "right",
                            align = "hv", 
                            font.label = list(size = 10, color = "black", family = NULL, position = "top"))
}

figure1 <- annotate_figure(figures[[1]],
                           top = text_grob("Linear regression: squared bias vs. k", size = 10, x = 0.4),
                           bottom = text_grob("k", size = 10, x = 0.4),
                           left = text_grob("Squared bias", size = 10, rot = 90))
figure2 <- annotate_figure(figures[[2]],
                           top = text_grob("Linear regression: variance vs. k", size = 10, x = 0.4),
                           bottom = text_grob("k", size = 10, x = 0.4),
                           left = text_grob("Variance", size = 10, rot = 90))
figure3 <- annotate_figure(figures[[3]],
                           top = text_grob("Linear regression: mean squared error (MSE) vs. k", size = 10, x = 0.4),
                           bottom = text_grob("k", size = 10, x = 0.4),
                           left = text_grob("MSE", size = 10, rot = 90))

counts <- list()
corrects <- c()
optplotsout <- list()
pred <- list()
library(tidyr)
for(i in 1:length(n)){
  
  # read results
  theta.hat <- mse.hat <- array(unlist(out[[i]]), dim = c(length(k[[i]]), n_mod, nsim))
  results <- matrix(nrow = nsim, ncol = length(k[[i]]))
  for(sim in 1:nsim){
    results[sim, ] <- apply(theta.hat[,, sim], 1, which.min)
  }
  corrects <- counts[[i]] <- colSums(results == 2)
  
  # fit
  lm.fit <- lm(corrects ~ k + I(1/k), data = data.frame(corrects = corrects, k = k[[i]]))
  lm.pred <- predict(lm.fit, data.frame(k = 2:n[[i]]))
  pred[[i]] <- as.numeric(lm.pred)
  
  # plot
  optplotsout[[i]] <- ggplot(NULL)+
    geom_point(data = data.frame(f = k[[i]], corrects), aes(x = f, y = corrects), size = 0.1) +
    geom_line(data = data.frame(f = 2:n[[i]], lm.pred), aes(x = f, y = lm.pred), size = 0.4, color = "black") +
    ggtitle(paste(titles[i])) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(colour="black", size=0.4),
          axis.text=element_text(size=10, color="black"),
          plot.title=element_text(hjust=0.5,size=10, color = "black")) +
    geom_vline(xintercept = sqrt(n[[i]]), color = "red", size = 0.4)
}

optfigures <- ggarrange(optplotsout[[1]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[2]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[3]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[4]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[5]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[6]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[7]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[8]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[9]] + rremove("ylab") + rremove("xlab"),
                        optplotsout[[10]] + rremove("ylab") + rremove("xlab"),
                        labels = NULL,
                        ncol = 2, nrow = 5,
                        common.legend = TRUE,
                        legend = "right",
                        align = "hv", 
                        font.label = list(size = 10, color = "black", family = NULL, position = "top"))

figure7 <- annotate_figure(optfigures,
                           top = text_grob("Linear regression: lowest-error model selection", size = 10),
                           bottom = text_grob("k", size = 10),
                           left = text_grob("Number of times true model had lowest MSE", size = 10, rot = 90))

library("smerc")
ep.linear <- c()
for(i in 1:length(k)){
  ep.linear[i] <- elbow_point(2:n[[i]], pred[[i]])$x
}

samp.n.linear <- seq(100, 1000, 100)
opt.k.linear <- c(14, 20, 25, 28, 32, 35, 37, 40, 42, 45)
opt.k.prop.linear <- opt.k.linear/samp.n.linear

samp.n.lasso <- c(250, 500, 750, 1000)
opt.k.lasso <- c(22, 32, 39, 45)
opt.k.prop.lasso <- opt.k.lasso/samp.n.lasso

opt.k.raw <- ggplot(data.frame(f = samp.n.linear, opt.k.linear), aes(x = f, y = opt.k.linear))+
  geom_line(data = data.frame(f = min(n):max(n), g = sqrt(2*(min(n):max(n)))), aes(x = f, y = g), size = 0.4, color = "black") +
  geom_point(aes(color = "Linear regression"), size = 0.4) +
  geom_point(data = data.frame(f = samp.n.lasso, opt.k.lasso), aes(x = f, y = opt.k.lasso, color = "LASSO regression"), size = 0.4) +
  labs(x = NULL, y = NULL) +
  ggtitle(bquote(atop("Optimal fold number"~(k[optimal]^"*"), "vs. sample size (n)"))) +
  xlab("n") +
  ylab(bquote(k[optimal]^"*")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black", size=0.4),
        axis.text=element_text(size=10, color="black"),
        axis.title=element_text(size=10),
        plot.title=element_text(hjust=0.5,size=10, color = "black"),
        legend.position = c(1, 1),
        legend.justification = c("right", "top"),
        legend.key=element_rect(fill = NA),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black")) +
  scale_colour_manual(name = "Model", values = c("blue", "red"))

opt.k.prop <- ggplot(data.frame(f = samp.n.linear, opt.k.prop.linear), aes(x = f, y = opt.k.prop.linear))+
  geom_line(data = data.frame(f = min(n):max(n), g = sqrt(2*(min(n):max(n)))/(min(n):max(n))), aes(x = f, y = g), size = 0.4, color = "black") +
  geom_point(aes(color = "Linear regression"), size = 0.4) +
  geom_point(data = data.frame(f = samp.n.lasso, opt.k.prop.lasso), aes(x = f, y = opt.k.prop.lasso, color = "LASSO regression"), size = 0.4) +
  labs(x = NULL, y = NULL) +
  ggtitle(bquote(atop("Optimal fold number proportion of sample size"~(k[optimal]^"*"/n), "vs. sample size (n)"))) +
  xlab("n") +
  ylab(bquote(k[optimal]^"*"/n)) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black", size=0.4),
        axis.text=element_text(size=10, color="black"),
        axis.title=element_text(size=10),
        plot.title=element_text(hjust=0.5,size=10, color = "black"),
        legend.position = c(1, 1),
        legend.justification = c("right", "top"),
        legend.key=element_rect(fill = NA),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black")) +
  scale_colour_manual(name = "Model", values = c("blue", "red"))

figure9 <- ggarrange(opt.k.raw, opt.k.prop,
            labels = NULL,
            ncol = 1, nrow = 2,
            common.legend = TRUE,
            legend = "right",
            align = "hv", 
            font.label = list(size = 10, color = "black", family = NULL, position = "top"))
