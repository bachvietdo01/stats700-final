sample_post <- function(m, n_sample = 10000, iter = 1000, d = 2) {
  K = rep(0, iter)
  z = list()
  X = array(0, c(iter, n_sample, d))
  rhos = rep(0, iter)
  
  for(s in 1:iter) {
    post_samp <- m$generate_sample(n = n_sample)
    rhos[s] <- post_samp$rho
    K[s]  <- length(unique(post_samp$z))
    z[[s]] <- post_samp$z
    X[s,,] = post_samp$X
  }
  
  return(list(K = K, X = X, z = z, rho = rhos))
}

sample_theta <- function(k = 1, X, z, K) {
  d = dim(X[1,,])[2]
  iter = length(K)
  mu = array(0, c(sum(K == k), d, k))
  Sigma = array(0,c(sum(K == k), d, d, k))
  pi = matrix(0, nrow = sum(K == k), ncol = k)
  
  i = 1
  for(s in 1:iter) {
    if (K[s] == k) {
      
      assign <- z[[s]]
      for(t in unique(assign)) {
        Xt = X[s,assign == t,]
        mut = colMeans(matrix(Xt, ncol = 2))
        mu[i,,t] = mut
        Sigma[i,,,t] = cov(Xt)
        pi[i,t] = sum(assign == t) / length(assign)
      }
      i = i + 1
    }
  }
  
  return(list(pi = pi, mu = mu, Sigma = Sigma))
}

hist_rho <- function(rho, rho_true = NULL) {
  # Plot dth dim of pi
  par(mfrow = c(1, 1))
  hist(rho, xlab = "rho",  main="Signal ratio")
  if (!is.null(rho_true)) {
    abline(v = rho_true, col = "red", lwd = 2)
  }
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


set.seed(0)
post_samp = sample_post(chandra)
Kfreq = table(post_samp$K)
par(mfrow = c(1, 1))
hist(rep(6,1000), breaks = c(3,4,5,6,7,8), main = 'K clusters', xlab = 'K')
output <- sample_theta(k = 6, post_samp$X, post_samp$z, post_samp$K)
hist_rho(post_samp$rho)
hist_pi(output$pi)
hist_mu(output$mu, data$mu)
hist_Sigma(output$Sigma, data$Sigma)

# get credible intervable
quantile(output$pi[4], c(0.025, 0.975))



