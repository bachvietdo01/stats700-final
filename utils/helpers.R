logdet_chol = function(L) 2*sum(log(diag(L)))

loglik_marginal_NIW <- function(N, D, kappa, nu, S_chol, kappa0, nu0, S0_chol){
  -0.5*N*D*log(pi) - 0.5*D*(log(kappa) - log(kappa0)) - 0.5*(nu*logdet_chol(S_chol) - nu0*logdet_chol(S0_chol)) + sum(lgamma(0.5*(nu - 1:D + 1)) - lgamma(0.5*(nu0 - 1:D + 1)))
}

add_sample_and_get_cholS_helper <- function(x, m, kappa, nu, L){
  m2 <- m + (x - m) / (kappa + 1L)
  chol_update(L, sqrt((kappa + 1L)/kappa) * (x - m2))
}

softmax <- function(...){
  logx <- unlist(list(...))
  normx <- logx - max(logx)
  exp(normx) / sum(exp(normx))
}

get_KingComponent_example_init_pars = function(D) {
  return(list(D = D, kappa0 = 1, nu0 = D + 2,
              m0 = rep(0, D), Gamma0 = diag(D),
              eta0 = 1, eta = 1.5))
}

get_IMVNComponent_example_init_pars = function(D) {
  return(list(D = D, kappa0 = 1, nu0 = D + 2,
              m0 = rep(0, D), s0.2 = 1))
}

get_GMVNComponent_example_init_pars = function(D) {
  return(list(D = D, kappa0 = 1, nu0 = D + 2,
              m0 = rep(0, D), S0 = diag(D)))
}

get_GMVNComponent_example_init_pars = function(D) {
  kappa0 = 1
  nu0 = D + 2
  m0 <- rep(0, D)
  S0 <- diag(D) 
  L0 <- chol(S0 + kappa0 * outer(m0, m0))
  
  return(list(D = D, kappa0 = kappa0, nu0 = nu0,
              m0 = m0, S0 = S0, L0 = L0))
}

get_EMNormComponent_example_init_pars = function(D) {
  return(list(mu = rep(0, D), Sigma = diag(D)))
}

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
