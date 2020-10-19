require(data.table)

plot_2D_GMM = function(obj, z = NULL){
  X = obj$X
  D = obj$D
  
  if(is.null(z)){
    z = factor(obj$z)
  } else{
    z = factor(z)
  }
  
  # construct joint mu matrix
  mus = matrix(0, nrow = D, ncol = obj$K)
  for(k in 1: obj$K) {
    mus[,k] = obj$sig_components[[k]]$mu
  }
  
  # plot out the inference result
  df = data.frame(x = X[,1], y = X[,2], z = z)
  p <-   ggplot(data = df) + geom_point(aes(x,y, color=z)) +
    geom_point(data = data.frame(x = mus[1, ], y = mus[2, ]),
               aes(x,y), color = 'red') +
    theme_bw() + theme(legend.position = "none") + coord_fixed()
   
  for(k in 1:obj$K) {
    # get fitted density countour plot
    out = get_grid_fit_density(X, obj$w[k], obj$sig_components[[k]])
    grid = out$grid
    fd = out$fd
    
    # add to plot
    p = p + stat_contour(data =data.frame(grid, density = fd),
                 aes(x,y, z = density), bins = 4, color = 'blueviolet')
  }
  
  p
}

get_grid_fit_density <- function(X, pik, clusterk) {
  xmax = apply(X, 2, max)
  xmin = apply(X, 2, min)
  
  # create a 2D grind points
  grid = as.data.table(expand.grid(
    x = seq(from = xmin[1], to = xmax[1], length.out = 100),
    y = seq(from = xmin[2], to = xmax[2], length.out = 100)))
  
  # compute fitted gmm density 
  muk = clusterk$mu
  Sigmak = clusterk$Sigma
  fit_dens = exp(log(pik) +  mvtnorm::dmvnorm(grid, mean = muk, sigma = Sigmak, 
                                              log = T))
  
  return(list(grid = grid, fd = fit_dens))
}
