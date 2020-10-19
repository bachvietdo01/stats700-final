library(ggplot2)
setwd('~/Downloads/MFMproj/uniform_background/')

Rcpp::sourceCpp('./helpers.cpp')
source('_component.R')
source('_mfmmixture.R')
source('helpers.R')
source('visualisation.R')

# plot data
ggplot() + geom_point(data = data.frame(x = data$X[,1], y = data$X[, 2]), 
                      aes(x, y))

# prepare to run Gibbs
X = as.matrix(data$X)
K = 20

s = rep(1, nrow(X))
z = kmeans(x = X, K)$cluster
s = em$s
z = em$z


m0_init = matrix(0, ncol = 6, nrow = 2)
m0_init[, 1] = c(4045, 4158)
m0_init[, 2] = c(4045, 4185)
m0_init[, 3] = c(4055, 4140)
m0_init[, 4] = c(4055, 4150)
m0_init[, 5] = c(4055, 4175)
m0_init[, 6] = c(4070, 4180)


chandra = MFMMixture(K = K, D = 2, X = X, s = s, z = z, m0 = colMeans(X),
                 m0_init = m0_init)

# run MFM Gibbs
niters = 10
K_iters = rep(0, niters)
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
    print(paste(c('iter:', i)))
    print(paste(c('Number of clusters', chandra$K)))
    print(plot_2D_GMM(chandra))
  }
  
  chandra$collapsed_gibbs()
  K_iters[i] = chandra$K
  #chandra$merge_split()
  #chandra$tidy_up_clusters()
}