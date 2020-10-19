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
