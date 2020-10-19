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
