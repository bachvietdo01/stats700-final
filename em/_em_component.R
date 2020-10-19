require(mvtnorm)

EMNormComponent <- setRefClass("EMNormComponent",
                               fields = list(
                                 D = "integer",
                                 Sigma = "matrix",
                                 mu = "numeric"),
                               
                               methods = list(
                                 initialize = function(init_pars){
                                   mu <<- init_pars$mu
                                   Sigma <<- init_pars$Sigma
                                   D <<- length(mu)
                                 },
                                 M_step_update = function(X, z_hat) {
                                   N <- sum(z_hat)
                                   
                                   mu <<- colSums(X  * z_hat) / N
                                     
                                   X_bar <- sweep(X, 2, mu, '-') * sqrt(z_hat)
                                   if(sum(X_bar) != 0) {
                                     Sigma <<- t(X_bar) %*% X_bar / N
                                   }
                                    
                                   if(any(is.na(mu))) {
                                     stop("update mu")
                                   }
                                   if(any(is.na(Sigma)) || sum(Sigma) == 0) {
                                     stop("update Sigma")
                                   }
                                 },
                                 get_loglik = function(x) {
                                   return(mvtnorm::dmvnorm(x, mean = mu, sigma = Sigma, log = TRUE))
                                 },
                                 generate_sample = function(n) {
                                   return(mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma))
                                 }
                               ))
