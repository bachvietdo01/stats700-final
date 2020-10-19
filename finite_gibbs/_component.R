require(mvtnorm)

GMVNComponent <- setRefClass("GMVNComponent",
                               fields = list(
                                 N = "integer",
                                 D = "integer",
                                 kappa = "numeric",
                                 nu = "numeric",
                                 m = "numeric",
                                 L = "matrix",
                                 S = "matrix",
                                 kappa0 = "numeric",
                                 nu0 = "numeric",
                                 m0 = "numeric",
                                 L0 = "matrix",
                                 Sigma = "matrix",
                                 mu = "numeric"),
                               
                               methods = list(
                                 initialize = function(init_pars, X = NULL){
                                   D <<- as.integer(init_pars$D)
                                   kappa0 <<- init_pars$kappa0
                                   nu0 <<- init_pars$nu0
                                   m0 <<- init_pars$m0
                                   S0 <- init_pars$S0
                                   L0 <<- chol(S0 + kappa0 * outer(m0, m0))
                                   N <<- 0L
                                   kappa <<- kappa0
                                   nu <<- nu0
                                   m <<- m0
                                   L <<- L0
                                   S <<- diag(D) # ???: wrong
                                   
                                   if(!is.null(X)){
                                     N_k <- nrow(X)
                                     kappa_k <- kappa + N_k
                                     m_k <- (kappa*m + colSums(X)) / kappa_k
                                     S_k <- S + kappa * outer(m, m) + t(X) %*% X
                                     L <<- chol(S_k - kappa_k * outer(m_k, m_k))
                                     if(any(is.nan(L))) stop("init problem")
                                     N <<- N + N_k
                                     nu <<- nu + N_k
                                     kappa <<- kappa_k
                                     m <<- m_k
                                     S <<- S_k
                                   }
                                   
                                   update_IW_pars()
                                 },
                                 is_empty = function(){
                                   ifelse(N == 0, TRUE, FALSE)
                                 },
                                 add_sample = function(x){
                                   kappa <<- kappa + 1L
                                   m <<- ((kappa - 1)*m + x) / kappa
                                   nu <<- nu + 1L
                                   L <<- chol_update(L, sqrt(kappa/(kappa - 1)) * (x - m))
                                   if(any(is.nan(L))) stop("Add sample")
                                   N <<- N + 1L
                                   S <<- L %*% t(L)
                                   #S <<- S + outer(x, x)
                                   #testthat::expect_equal(chol(S - kappa * outer(m, m)), L)
                                   update_IW_pars()
                                 },
                                 rm_sample = function(x){
                                   L <<- chol_downdate(L, sqrt(kappa/(kappa - 1)) * (x - m))
                                   if(any(is.nan(L))){
                                     
                                     stop("Rm sample")
                                   }
                                   kappa <<- kappa - 1L
                                   m <<- ((kappa + 1)*m - x) / kappa
                                   nu <<- nu - 1L
                                   N <<- N - 1L
                                   S <<- L %*% t(L)
                                   #S <<- S - outer(x, x)
                                   #testthat::expect_equal(chol(S - kappa * outer(m, m)), L)
                                   # if(mean(abs(chol(S - kappa * outer(m, m)) - L)) > 1e-5) stop("We have a problem")
                                   update_IW_pars()
                                 },
                                 get_cholS = function(){
                                   # subtracts the cluster means from the covariance matrix
                                   # and returns Cholesky decomposition
                                   # LL <- chol(S - kappa * outer(m, m))
                                   return(L)
                                 },
                                 get_S = function(){
                                   #S - kappa * outer(m, m)
                                   return(L %*% t(L))
                                 },
                                 add_sample_and_get_cholS = function(x){
                                   # call the generic function
                                   return(add_sample_and_get_cholS_helper(x, m, kappa, nu, L))
                                 },
                                 marginal_loglik = function(){
                                   loglik_marginal_NIW(N, D, kappa, nu, L, kappa0, nu0, L0)
                                 },
                                 posterior_predictive = function(x){
                                   # prediction when data point x were in the current cluster
                                   L_updated <- add_sample_and_get_cholS(x)
                                   # calculate loglikelihood
                                   loglik <- loglik_marginal_NIW(1, D, kappa + 1L, nu + 1L, L_updated, kappa, nu, L)
                                   # loglik <- -0.5*D*log(pi) - 0.5*D*log((kappa+1)/kappa) - 0.5*(nu+1)*logdet_chol(L_updated) + lgamma(0.5*(nu + 1)) + 0.5*nu*logdet_chol(L_current) - lgamma(0.5*(nu + 1 - D))
                                   loglik
                                 },
                                 
                                 get_loglik = function(x) {
                                   return(mvtnorm::dmvnorm(x, mean = mu, sigma = Sigma, log = TRUE))
                                 },
                                 
                                 update_IW_pars = function(){
                                   #Gamma <- solve(get_S())
                                   #Sigma <<- solve(rWishart(1, nu, Gamma)[, , 1])
                                   Sigma <<- MCMCpack::riwish(nu, get_S())
                                   mu <<- as.numeric(mvtnorm::rmvnorm(1, mean = m, sigma = 1/kappa * Sigma))
                                 }
                               ))

IMVNComponent <- setRefClass("IMVNComponent",
                               fields = list(
                                 N = "integer",
                                 D = "integer",
                                 kappa = "numeric",
                                 nu = "numeric",
                                 m = "numeric",
                                 s.2 = "numeric",
                                 kappa0 = "numeric",
                                 nu0 = "numeric",
                                 m0 = "numeric",
                                 s0.2 = "numeric",
                                 mu = "numeric",
                                 tau = "numeric"),
                               
                               methods = list(
                                 initialize = function(init_pars, X = NULL){
                                   
                                   D <<- as.integer(init_pars$D)
                                   kappa0 <<- init_pars$kappa0
                                   nu0 <<- init_pars$nu0
                                   m0 <<- init_pars$m0
                                   s0.2 <<- init_pars$s0.2
                                   
                                   N <<- 0L
                                   kappa <<- kappa0
                                   nu <<- nu0
                                   m <<- m0
                                   s.2 <<- s0.2
                                   
                                   if(!is.null(X) && dim(X)[1] != 0){
                                     N_k <- nrow(X)
                                     kappa_k <- kappa + N_k
                                     m_k <- (kappa*m + colSums(X)) / kappa_k
                                     
                                     if(N_k > 1) {
                                       xbar <- colMeans(X)
                                     } else {
                                       xbar <- X
                                     }
                                     
                                     Xbar = sweep(X, 2, xbar)
                                     s_k.2 <- nu * s.2 + sum(diag(Xbar %*% t(Xbar))) + N_k * kappa /(kappa_k) * t(m - xbar) %*% (m - xbar)
                                     s_k.2 <- s_k.2 / (nu + D * N_k)

                                     N <<- N + N_k
                                     nu <<- nu + D * N_k
                                     kappa <<- kappa_k
                                     m <<- m_k
                                     s.2 <<- as.numeric(s_k.2)
                                   }
                                   
                                   update_pars()
                                 },
                                 is_empty = function(){
                                   ifelse(N == 0, TRUE, FALSE)
                                 },
                                 add_sample = function(x){
                                   kappa <<- kappa + 1L
                                   nu <<- nu + D
                                   s.2 <<- as.numeric((nu -D) * s.2 + (kappa - 1) / kappa * t(m - x) %*% (m - x))
                                   s.2 <<- s.2 / nu
                                   m <<- ((kappa - 1)*m + x) / kappa
                                   
                                   N <<- N + 1L
                                   update_pars()
                                 },
                                 rm_sample = function(x){
                                   kappa <<- kappa - 1L
                                   nu <<- nu - D
                                   m <<- ((kappa + 1)*m - x) / kappa
                                   s.2 <<- as.numeric((nu + D) * s.2 - kappa/ (kappa + 1) * t(m -x) %*% (m -x))
                                   s.2 <<- s.2 / nu
                                   
                                   N <<- N - 1L
                                   update_pars()
                                 },
                                 
                                 get_loglik = function(x) {
                                   return(mvtnorm::dmvnorm(x, mean = mu, sigma = tau * diag(D), log = TRUE))
                                 },
                                 
                                 update_pars = function(){
                                   tau <<- (nu * s.2)/rchisq(1, df = nu)
                                   mu <<- as.numeric(mvtnorm::rmvnorm(1, mean = m, sigma = tau /kappa * diag(D)))
                                 }
                               ))

KingComponent <- setRefClass("KingComponent",
                             fields = list(
                               N = "integer",
                               D = "integer",
                               mu = "numeric",
                               S.L = "matrix",
                               Sigma.L = "matrix",
                               Psi = "matrix",
                               Sigma = "matrix",
                               eta = "numeric",
                               m = "numeric",
                               kappa = "numeric",
                               nu = "numeric",
                               S = "matrix",
                               S.L0 = "matrix",
                               m0 = "numeric",
                               kappa0 = "numeric",
                               nu0 = "numeric",
                               eta0 = "numeric",
                               Gamma0 = "matrix"),
                             methods = list(
                               initialize = function(init_pars, X = NULL){
                                 D <<- as.integer(init_pars$D)
                                 
                                 eta <<- init_pars$eta
                                 eta0 <<- init_pars$eta0
                                 Gamma0 <<- init_pars$Gamma0
                                 
                                 m0 <<- init_pars$m0
                                 kappa0 <<- init_pars$kappa0
                                 nu0 <<- 2 * init_pars$eta - 1
                                 
                                 Sigma.L <<- CholWishart::rCholWishart(1, nu0, Gamma0)[,,1]
                                 Sigma <<- t(Sigma.L) %*% Sigma.L
                                 
                                 S <<-Sigma
                                 S.L0 <<- chol(S + kappa0 * outer(m0, m0))
                                 
                                 N <<- 0L
                                 kappa <<- kappa0
                                 nu <<- nu0
                                 m <<- m0
                                 
                                 if(!is.null(X)){
                                   N_k <- nrow(X)
                                   kappa_k <- kappa + N_k
                                   m_k <- (kappa*m + colSums(X)) / kappa_k
                                   S_k <- S + kappa * outer(m, m) + t(X) %*% X
                                   S.L <<- chol(S_k - kappa_k * outer(m_k, m_k))
                                   
                                   if(any(is.nan(S.L))) stop("init problem")
                                   
                                   N <<- N + N_k
                                   nu <<- nu + N_k
                                   kappa <<- kappa_k
                                   m <<- m_k
                                   S <<- S_k
                                 }
                                 
                                 update_King_pars()
                               },
                               is_empty = function(){
                                 ifelse(N == 0, TRUE, FALSE)
                               },
                               add_sample = function(Xk){
                                 block_update(Xk)
                               },
                               rm_sample = function(Xk){
                                 block_update(Xk)
                               },
                               block_update = function(Xk) {
                                 N <<- 0L
                                 kappa <<- kappa0
                                 nu <<- nu0
                                 m <<- m0
                                 
                                 S <<- Sigma
                                
                                 N_k <- nrow(Xk)
                                 kappa_k <- kappa + N_k
                                 m_k <- (kappa*m + colSums(Xk)) / kappa_k
                                 S_k <- S + kappa * outer(m, m) + t(Xk) %*% Xk
                                 S.L <<- chol(S_k - kappa_k * outer(m_k, m_k))
                                 
                                 if(any(is.nan(S.L))) stop("init problem")
                                 
                                 N <<- N + N_k
                                 nu <<- nu + N_k
                                 kappa <<- kappa_k
                                 m <<- m_k
                                 S <<- S_k
                                 
                                 update_King_pars()
                               },
                               # add_sample = function(x){
                               #   kappa <<- kappa + 1L
                               #   m <<- ((kappa - 1)*m + x) / kappa
                               #   
                               #   nu <<- nu + 1L
                               #   S.L <<- chol_update(S.L, sqrt(kappa/(kappa - 1)) * (x - m))
                               #   
                               #   if(any(is.nan(S.L))) stop("Add sample")
                               #   
                               #   N <<- N + 1L
                               #   S <<- t(S.L) %*% S.L
                               #   
                               #   Sigma_old <- Sigma
                               #   # get new mu, Sigma, Psi
                               #   update_King_pars()
                               #   
                               #   # update S for new iteration
                               #   S <<- S - Sigma_old + Sigma
                               #   S.L <<- chol(S)
                               #   if(any(is.nan(S.L))){
                               #     stop("add sample")
                               #   }
                               # },
                               # rm_sample = function(x){
                               #   S.L <<- chol_downdate(S.L, sqrt(kappa/(kappa - 1)) * (x - m))
                               #   if(any(is.nan(S.L))){
                               #     stop("Rm sample")
                               #   }
                               #   kappa <<- kappa - 1L
                               #   m <<- ((kappa + 1)*m - x) / kappa
                               #   nu <<- nu - 1L
                               #   N <<- N - 1L
                               #   S <<- t(S.L) %*% S.L
                               #   #S <<- S - outer(x, x)
                               #   #testthat::expect_equal(chol(S - kappa * outer(m, m)), L)
                               #   # if(mean(abs(chol(S - kappa * outer(m, m)) - L)) > 1e-5) stop("We have a problem")
                               #   
                               #   Sigma_old <- Sigma
                               #   # get new mu, Sigma, Psi
                               #   update_King_pars()
                               #   S <<- S - Sigma_old + Sigma
                               #   S.L <<- chol(S)
                               #   if(any(is.nan(S.L))){
                               #     stop("rm sample")
                               #   }
                               # },
                               get_cholS = function(){
                                 # subtracts the cluster means from the covariance matrix
                                 # and returns Cholesky decomposition
                                 # LL <- chol(S - kappa * outer(m, m))
                                 return(S.L)
                               },
                               get_S = function(){
                                 #S - kappa * outer(m, m)
                                 return(t(S.L) %*% S.L)
                               },
                               get_loglik = function(x) {
                                 return(mvtnorm::dmvnorm(x, mean = mu, sigma = Psi, log = TRUE))
                               },
                               update_King_pars = function(){
                                 Psi <<- CholWishart::rInvWishart(1, nu, S)[,,1]
                                 if(any(is.nan(Psi))){
                                   stop("Stop update king pars")
                                 }
                                 mu <<- as.numeric(mvtnorm::rmvnorm(1, mean = m, sigma = 1/kappa * Psi))
                                 mu <<- m0
                                 
                                 # update Sigma
                                 Gamma <- Gamma0 - Gamma0 %*% solve(Psi + Gamma0) %*% Gamma0
                                 Sigma.L <<- CholWishart::rCholWishart(1, 2 * eta - 1 + eta0, Gamma)[,,1]
                                 if(any(is.nan(Sigma.L))){
                                   stop("Stop update king pars")
                                 }
                                 Sigma <<- t(Sigma.L) %*% Sigma.L
                               }
                             ))
