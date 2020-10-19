library(ggplot2)
setwd('~/Documents/MixSigBkg/R/em/')

Rcpp::sourceCpp('./helpers.cpp')
source('../visualisation.R')
source('_em_component.R')
source('_em_mixture.R')
source('_em_mmsignal.R')
source('helpers.R')

# plot data
data <- list(X = spatial[1:1000,])
ggplot() + geom_point(data = data.frame(x = data$X[,1], y = data$X[,2]), aes(x, y)) 

# prepare to run Gibbs
X = as.matrix(data$X)
Ks = 6
s = rep(1, nrow(X))
z = kmeans(x = X, Ks + 1)$cluster

m_init = matrix(rep(apply(X, 2, min), Ks), nrow = dim(X)[2] )
m_init[, 1] = c(4042, 4155)
m_init[, 2] = c(4045, 4180)
m_init[, 3] = c(4055, 4140)
m_init[, 4] = c(4055, 4150)
m_init[, 5] = c(4055, 4175)
m_init[, 6] = c(4070, 4180)


gem_init_pars = get_EMNormComponent_example_init_pars(2)
gem_init_pars$m_init = m_init

chandra = EMMixSigAndBg(Ks = Ks, D = 2, X = X, s = s, z = z,
               sig_componentClass = EMNormComponent, sig_init_pars = gem_init_pars)
  
# run EM
niters = 20

loglik <- rep(0, niters)
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
    print(paste(c('iter:', i)))
    print(plot_2D_MM_signal(chandra))
  }
    
  chandra$em_update()
  loglik[i] <- chandra$get_loglik()
}

plot_em_loglik(loglik)

s = chandra$signal[[1]]$s
z = chandra$signal[[1]]$z

rm(list=ls()[! ls() %in% c("data","s", "z")])
