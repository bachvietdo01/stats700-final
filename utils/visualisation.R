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

plot_result = function(X, z) {
  df = data.frame(x = X[,1], y = X[,2], z = as.factor(z))
  
  p <- ggplot(df) + geom_point(aes(x,y, col = z)) + theme_bw()
  p
}

hist_rho <- function(rho, rho_true = NULL) {
  # Plot dth dim of pi
  par(mfrow = c(1, 1))
  hist(rho, xlab = "rho",  main="Signal ratio")
  if (!is.null(rho_true)) {
    abline(v = rho_true, col = "red", lwd = 2)
  }
}

vb_hist_rho <- function(rho, M = 1, s = 1, rho.true = NULL, nbins = 30) {
  par(mfrow = c(1, 1))
  # Plot rho pdf
  rho.crop = rho[-(1:M)]
  
  xlim_lb = min(rho.crop)
  xlim_ub = max(rho.crop)
  
  if (!is.null(rho.true)) {
    xlim_lb = min(xlim_lb, rho.true)
    xlim_ub = max(xlim_ub, rho.true)
  }
  
  
  hist(rho.crop[seq(1, length(rho.crop), s)], nbins, xlab = 'rho', main =
         'Signal ratio', xlim = c(xlim_lb, xlim_ub))
  if (!is.null(rho.true)) {
    abline(v = rho.true, col = "red", lwd = 2)
  }
  effectiveSize(rho.crop[seq(1, length(rho.crop), s)])
}

hist_pi <- function(pi, pi_true = NULL) {
  K = dim(pi)[2]
  par(mfrow = c(1, K))
  # Plot dth dim of pi
  for (k in 1:K) {
    pi.crop = pi[,k]
    hist(pi.crop, 15, xlab = paste(c('pi[', k,']')), 
         main=paste(c('Cluster ',k, ' weight')))
    if (!is.null(pi_true)) {
      abline(v = pi_true[k], col = "red", lwd = 2)
    }
  }
}

hist_mu <- function(mu, mu_true = NULL) {
  K = dim(mu[1,,])[2]
  d = dim(mu[1,,])[1]
  
  par(mfrow = c(d, K))
  
  for (p in 1:d) {
    for (k in 1:K) {
      mu.cropik = mu[,p,k]
      
      xlim_lb = min(mu.cropik)
      xlim_ub = max(mu.cropik)
      
      if (!is.null(mu_true)) {
        xlim_lb = min(xlim_lb, mu_true[p,k])
        xlim_ub = max(xlim_ub, mu_true[p,k])
      }
      
      hist(mu.cropik, 15, xlab = paste(c('mu[', p, ',', k,']'),
                                       collapse = ''), main = 'mu[d, k]', xlim = c(xlim_lb,
                                                                                   xlim_ub))
      if (!is.null(mu_true)) {
        abline(v = mu_true[p, k], col = "red", lwd = 2)
      }
    }
  }
}

hist_Sigma <- function(Sigma, Sigma_true = NULL) {
  K = dim(Sigma[1,,,])[3]
  d = dim(Sigma[1,,,])[1]
  
  par(mfrow = c(K, d * d))
  
  for (k in 1:K) {
    for (p1 in 1:d) {
      for (p2 in 1:d) {
        Sigma.cropijk = Sigma[,p1, p2, k]
        
        hist(Sigma.cropijk, 15, xlab = paste(c('Sigma[', p1, ',', p2,',', k,']')
                                             ,collapse = ''),
             main = 'Sigma[d,d, k]')
        
        if (!is.null(Sigma_true)) {
          abline(v = Sigma_true[p1, p2, k], col = "red", lwd = 2)
        }
      }
    }
  }
}

make_visual_obj_from_vb <- function(vb_run, X) {
  K = dim(vb_run$r_nk)[2] - 1
  d = dim(X)[2]
  z <- apply(vb_run$r_nk, 1, function(p) sample(1:length(p), 1, prob = p))
  z[z == (K + 1)] <- 0
  s <- z
  s[s != 0] <- 1
  
  res = list()
  res$z = z
  res$s = s
  res$X = data.frame(X = X[,1], Y = X[,2])
  res$K = K
  res$bg_lower = vb_run$bg_lower
  res$bg_upper = vb_run$bg_upper
  
  Sigma = array(0, c(d, d, K, 1))
  mu = array(0, c(d, K, 1))
  
  for(k in 1:K) {
    mu[, k, 1] = vb_run$m[, k]
    Sigma[,, k, 1] = solve(vb_run$W[,, k] * vb_run$nu[k])
  }
  
  res$mu = mu
  res$Sigma = Sigma
  
  return(res)
}

visualize_result <- function(res, X, K, Umin, Umax, true_mu = NULL, true_z = NULL, true_s = NULL) {
  s = res$s
  z = res$z
  
  
  pc = c('lightskyblue', 'gold', 'chartreuse4', 'magenta', 'tan1','mediumorchid',
         'blue', 'darkorange', 'black', 'cyan', 'indianred', 'lightseagreen',
         'pink', 'salmon3' )
  
  data.bg = as.data.frame(X[s == 0, ])
  colnames(data.bg) <- c("x", "y")
  
  p <- ggplot() + coord_fixed() + theme_bw() + 
    geom_point(data = data.bg, aes(x, y), colour = "grey50")
  
  for(k in 1:K) {
    data.ck = as.data.frame(X[z == k, ])      
    colnames(data.ck) <- c("x", "y")
    p <- p + geom_point(data = data.ck, aes(x, y), colour = pc[k])
  }
  
  
  grid <- as.data.table(expand.grid(x = seq(from = Umin[1], to = Umax[1], length.out = 100),
                                    y = seq(from = Umin[2], to = Umax[2], length.out = 100)))
  
  probs = rep(0, K)
  last_iter = dim(res$mu)[3]
  
  for(k in 1:K) {
    p <- p + geom_point(data = data.frame(X = res$mu[1,k,last_iter],
                                          Y = res$mu[2,k,last_iter]),
                        aes(X,Y), colour = 'red')
    
    probs = mvtnorm::dmvnorm(grid, res$mu[, k, last_iter], res$Sigma[,, k,last_iter])
    contour_probs <- as.data.frame(cbind(grid, density = probs))
    p <- p + geom_contour(data = contour_probs, aes(x,y,z = density),
                          colour = pc[k], bins = 5)
  }
  
  p
}

vb_get_sample_rho <- function(vb_run, N = 1000) {
  rhos <- rbeta(n = N, shape1 = vb_run$a[1], shape2 = vb_run$a[2])
  
  return(rhos)
}

vb_get_sample_pi <- function(vb_run, N = 1000) {
  pis <- rdirichlet(N, vb_run$alpha)
  
  return(t(pis))
}


vb_get_sample_mu <- function(vb_run, N = 1000) {
  d = dim(vb_run$m)[1]
  K = dim(vb_run$m)[2]
  
  mu_sample <- array(0, c(d, K, N))
  
  for(k in 1:K) {
    Gamma_k = vb_run$W[,, k] * vb_run$nu[k]
    mu_sample[,k,] =t(rmvnorm(N, vb_run$m[, k], solve(vb_run$beta[k] * Gamma_k)))
  }
  
  return(mu_sample)
}


vb_get_sample_Sigma <- function(vb_run, N = 1000) {
  d = dim(vb_run$m)[1]
  K = dim(vb_run$m)[2]
  
  sigma_sample <- array(0, c(d, d, K, N))
  
  for(k in 1:K) {
    gamma_k_sample <- rWishart(N, vb_run$nu[k], vb_run$W[,,k])
    
    for(i in 1:N) {
      sigma_sample[,,k,i] = solve(gamma_k_sample[,,i])
    }
  }
  
  return(sigma_sample)
}

