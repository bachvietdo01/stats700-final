library(ggplot2)
setwd('~/Downloads/MFMproj/uniform_background/')

Rcpp::sourceCpp('./helpers.cpp')
source('_component.R')
source('_mfmmixture.R')
source('_dpmixture.R')
source('helpers.R')
source('visualisation.R')

# plot data
indices = sample(1:nrow(data), 1500)
X <- data[indices,]
ggplot() + geom_point(data = data.frame(x = X[indices,1], y = X[indices, 2]), 
                      aes(x, y))

# prepare to run Gibbs
X = as.matrix(X)
K = 10
s = rep(1, nrow(X))
z = kmeans(x = X, K)$cluster

xmm = MFMMixture(K = K, D = 2, X = X, s = s, z = z, m0 = colMeans(X))

# run MFM Gibbs
niters = 50
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
    print(paste(c('iter:', i)))
    print(paste(c('Number of clusters', xmm$K)))
    print(plot_2D_GMM(xmm))
  }
  
  xmm$collapsed_gibbs()
  xmm$merge_split()
  xmm$tidy_up_clusters()
}