require(mvtnorm)

VBNormComponent <- setRefClass("VBNormComponent",
                               fields = list(
                                 D = "integer",
                                 Sigma = "matrix",
                                 mu = "numeric",
                                 W0 = "matrix",
                                 W0_inv = "matrix",
                                 m0 = "numeric",
                                 kappa0 = "numeric",
                                 nu0 = "numeric",
                                 W = "matrix",
                                 m = "numeric",
                                 kappa = "numeric",
                                 nu = "numeric",
                                 log_Lambda = "numeric",
                                 log_lb = "numeric"),
                               methods = list(
                                 initialize = function(init_pars){
                                   D <<- length(m0)
                                   
                                   # store inital values
                                   m0 <<- init_pars$m0
                                   W0 <<- init_pars$W0
                                   kappa0 <<- init_pars$kappa0
                                   nu0 <<- init_pars$nu0
                                   
                                   W0_inv <<- solve(W0)
                                   
                                   # assign current values
                                   m <<- m0
                                   W <<- W0
                                   kappa <- kappa0
                                   nu <<- nu0
                                   log_Lambda <<- sum(digamma((nu + 1 - c(1:D))/2)) +
                                     D*log(2) + log(det(W))
                                 },
                                 vb_Estep_update = function(X) {
                                   Xbar <- sweep(X,  2,  m, "-")
                                   log_q_nk <- D/2 *log(2*pi) + 0.5*log_Lambda -
                                     0.5*(D/kappa) - 0.5*nu * diag(Xbar %*% W %*% t(Xbar)) # log of 10.67
                                   
                                   return(log_q_nk)
                                 },
                                 vb_Mstep_update = function(X, r_nk) {
                                   Nk <- sum(r_nk) + 1e-10  # 10.51
                                    
                                   x_bar <- r_nk %*% X / Nk   # 10.52
                                   x_cen <- sweep(X, 2, x_bar, "-")
                                   S   <- t(x_cen) %*% (x_cen * r_nk) / Nk  # 10.53
                                   
                                   kappa <<- kappa0 + Nk   # 10.60
                                   nu <<- nu0 + Nk  # 10.63
                                   m <<- (1/kappa) * (kappa0*m0 + Nk * x_bar) # 10.61
                                   W <<- W0_inv + Nk * S +
                                       ((kappa0*Nk)/(kappa0 + Nk)) *tcrossprod((x_bar - m0)) # 10.62
                                   W <<- solve(W)
                                   
                                   log_Lambda <<- sum(digamma((nu + 1 - 1:D)/2)) +
                                       D*log(2) + log(det(W)) # 10.65
                                   
                                   # compute lower bound of p(x | rest)
                                   log_rho = digamma(a) - digamma(sum(a))
                                   log_pi <- digamma(alpha) - digamma(sum(alpha))
                                   
                                   lb_px <- Nk * (log_Lambda - D/kappa - nu *
                                                    matrix.trace(S %*% W) -
                                                    nu *t(x_bar - m) %*% W %*% (x_bar - m) -
                                                    D*log(2*pi)) # 10.71                                   
                                   lb_pml <- D*log(kappa0/(2*pi)) + log_Lambda -
                                     (D*kappa0)/kappa - kappa0*nu*t(m - m0) %*% W %*% (m - m0) 
                                   lb_pml2 <- nu * matrix.trace(W0_inv %*% W)
                                   lb_pml <- 0.5*lb_pml - 0.5*lb_pml2 + logB(W = W0,nu = nu0) + 
                                             0.5*(nu0 - D - 1) * log_Lambda # 10.74
                                   lb_qml <- 0.5*log_Lambda + 0.5*D*log(kappa/(2*pi)) - 0.5*D - 
                                     (-logB(W = W, nu = nu) - 0.5*(nu - D - 1)*log_Lambda + 0.5*nu*D) # 10.77
                                   

                                   log_lb <<- 0.5 * lb_px + lb_pml - lb_qml
                                   
                                   update_pars()
                                 },
                                 update_pars = function() {
                                   # draw the precision matrix
                                   L_Gamma <- CholWishart::rCholWishart(1, df = nu, Sigma = W)
                                   Gamma <- t(L_Gamma) %*% L_Gamma
                                   
                                   L_Sigma <- chol2inv(L_Gamma)
                                   Sigma <<- t(L_Sigma) %*% L_Sigma
                                   m <<- mvtnorm::rmvnorm(1, mean = m, sigma = Sigma)
                                 },
                                 get_log_lb = function() {
                                   return(log_lb)
                                 }
                                 get_loglik = function(x) {
                                   return(mvtnorm::dmvnorm(x, mean = mu, sigma = Sigma, log = TRUE))
                                 },
                                 generate_sample = function(n) {
                                   return(mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma))
                                 }
                               ))
