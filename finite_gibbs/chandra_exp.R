library(ggplot2)
setwd('~/Documents/MixSigBkg/R/')

Rcpp::sourceCpp('./helpers.cpp')
source('_component.R')
source('_mixture.R')
source('_fmmsignal.R')
source('helpers.R')
source('visualisation.R')

# plot data
data <- list(X = spatial[1:1000,])
ggplot() + geom_point(data = data.frame(x = data$X[,1], y = data$X[,2]), aes(x, y)) 

# prepare to run Gibbs
set.seed(0)
X = as.matrix(data$X)
Ks = 6
Kb = 1
#s = rep(1, nrow(X))

sz = kmeans(x = X, Ks)$cluster
bz = rep(0, nrow(X))

gmm_init_pars = get_GMVNComponent_example_init_pars(2)
m_init = matrix(rep(apply(X, 2, min), Ks), nrow = dim(X)[2] )
m_init[, 1] = c(4045, 4158)
m_init[, 2] = c(4045, 4185)
m_init[, 3] = c(4055, 4140)
m_init[, 4] = c(4055, 4150)
m_init[, 5] = c(4055, 4175)
m_init[, 6] = c(4070, 4180)

gmm_init_pars$m_init = m_init
gmm_init_pars$S0 = 9 * gmm_init_pars$S0
  
# initialize with em run
s = s
sz = z

chandra = Mixture(Ks = Ks, Kb = Kb, D = 2, X = X, s = s, sz = sz, bz = bz,
               sig_componentClass = GMVNComponent, sig_init_pars = gmm_init_pars)
  
# run Gibbs
niters = 5
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
      print(paste(c('iter:', i)))
      print(plot_2D_GMM_signal(chandra, chandra$signal[[1]]$z))
      for(k in 1:Ks) {
        print(paste(c('Signal cluster k size:', chandra$signal[[1]]$components[[k]]$N)))
      }
    }
    
  chandra$gibbs()
}