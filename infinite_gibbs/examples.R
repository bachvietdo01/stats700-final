library(ggplot2)
setwd('~/Downloads/MFMproj/uniform_background/')

attach('astronomydata/XMM - 2 sources/subset and plot/xmm_subset.RData',
                  name = 'xmm_data')
rm(list=setdiff(ls(), "data"))
data = spatial
data = data
detach('xmm_data')

Rcpp::sourceCpp('./helpers.cpp')
source('_component.R')
source('_mfmmixture.R')
source('_dpmixture.R')
source('helpers.R')
source('visualisation.R')

# plot data
indices = sample(1:nrow(data$X), 1000)
ggplot() + geom_point(data = data.frame(x = data$X[indices,1], y = data$X[indices, 2]), 
                      aes(x, y))

# prepare to run Gibbs
X = as.matrix(data$X)
K = 10
s = rep(1, nrow(X))
z = kmeans(x = X, K)$cluster
  
# create mixture models
m0_init = matrix(rep(apply(X, 2, min), K), nrow = dim(X)[2] )
m0_init[, 1] = c(4045, 4180)
m0_init[, 2] = c(4060, 4075)
m0_init[, 3] = c(4070, 4180)
m0_init[, 4] = c(4045, 4155)
m0_init[, 5] = c(4050, 4060)
m0_init[, 6] = c(4055, 4140)
  
m1 = MFMMixture(K = K, D = 2, X = X, s = s, z = z, m0 = colMeans(X))
  
# run MFM Gibbs
niters = 50
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
    print(paste(c('iter:', i)))
    print(paste(c('Number of clusters', m1$K)))
    print(plot_2D_GMM(m1))
  }
    
  m1$collapsed_gibbs()
  m1$merge_split()
  m1$tidy_up_clusters()
}


# Run DP Gibbs
X = as.matrix(data$X)
K = 10
s = rep(1, nrow(X))
z = kmeans(x = X, K)$cluster

m2 = DPMixture(K = K, D = 2, X = X, s = s, z = z, m0 = colMeans(X))

# run Gibbs
niters = 50
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
    print(paste(c('iter:', i)))
    print(paste(c('Number of clusters', m2$K)))
    print(plot_2D_GMM(m2))
  }
  
  m2$collapsed_gibbs()
  m2$merge_split()
}
