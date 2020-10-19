require(data.table)

plot_2D_GMM_signal = function(obj, z = NULL, scale_dens = 1.0){
  X = obj$X
  D = obj$D
  
  if(is.null(z)){
    z = factor(obj$z)
  } else{
    z = factor(z)
  }
  
  signal = obj$signal[[1]]
  
  # construct joint mu matrix
  mus = matrix(0, nrow = D, ncol = signal$K)
  for(k in 1: obj$Ks) {
    mus[,k] = signal$components[[k]]$mu
  }
  
  # plot out the inference result
  df = data.frame(x = X[,1], y = X[,2], z = z)
  p <-   ggplot(data = df) + geom_point(aes(x,y, color=z)) +
    geom_point(data = data.frame(x = mus[1, ], y = mus[2, ]),
               aes(x,y), color = 'red') +
    theme_bw() + theme(legend.position = "none")
  
  for(k in 1:obj$Ks) {
    # get fitted density countour plot
    out = get_grid_fit_density(X, signal$w[k], signal$components[[k]])
    grid = out$grid
    fd = scale_dens * out$fd
    
    # add to plot
    p = p + stat_contour(data =data.frame(grid, density = fd),
                         aes(x,y, z = density), bins = 4, color = 'royalblue')
  }
  
  p + coord_fixed(ratio = 1.0)
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

plot_2D_MM_signal = function(obj, z = NULL){
  X = obj$X
  if(is.null(z)){
    z = factor(obj$signal[[1]]$z)
  } else{
    z = factor(z)
  }
  df = data.frame(x = X[, 1], y = X[, 2], z = z)
  ggplot(df, aes(x, y, col = z)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "none")
}

plot_2D_MM_background = function(obj, z = NULL){
  X = obj$X
  if(is.null(z)){
    z = factor(obj$background[[1]]$z)
  } else{
    z = factor(z)
  }
  df = data.frame(X1 = X[, 1], X2 = X[, 2], z = z)
  ggplot(df, aes(X1, X2, col = z)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "none")
}

plot_em_loglik = function(loglik) {
  ggplot(data = data.frame(iters = 1: length(loglik), loglik = loglik)) + 
    geom_point(aes(iters,loglik)) + geom_line(aes(iters, loglik))
}

