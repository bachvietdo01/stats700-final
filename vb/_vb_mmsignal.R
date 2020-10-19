VBFMMSignal <- setRefClass("VBFMMSignal",
                       fields = list(
                         components = "list",
                         K = "integer",
                         D = "integer",
                         X = "matrix",
                         r_nk = "matrix",
                         log_r_nk = "matrix",
                         log_q_nk = "matrix",
                         ComponentClass ="refObjectGenerator",
                         init_pars = "list",
                         N = "integer",
                         bg_density = "numeric",
                         rho = "numeric",
                         w = "numeric",
                         a0 = "numeric",
                         alpha0 = "numeric",
                         a = "numeric",
                         alpha = "numeric",
                         log_w = "numeric",
                         log_rho = "numeric",
                         s = "numeric",
                         z = "numeric"),

                       methods = list(
                         initialize = function(K, D, X, s, z, 
                                               alpha0 = 1 / K,
                                               a0 = 0.5,
                                               bg_density,
                                               ComponentClass = NULL,
                                               init_pars = NULL) {
                           
                           testthat::expect_true(!is.null(ComponentClass))
                           
                           ComponentClass <<- ComponentClass
                           
                           D <<- as.integer(D)
                           K <<- as.integer(K)
                           X <<- X
                           
                           N <<- nrow(X)
                           s <<- s
                           z <<- z
                           
                           r_nk <<- matrix(0, nrow = N, ncol = K + 1)
                           log_r_nk <<- matrix(0, nrow = N, ncol = K + 1)
                           log_q_nk <<- matrix(0, nrow = N, ncol = K + 1)
                           
                           bg_density <<- bg_density
                           m0_inits <- kmeans(x = X, K)$centers
                           
                           for(k in 1:K){
                             inik <- init_pars
                             
                             if(!is.null(init_pars$m_init)) {
                               inik$mu <- init_pars$m_init[,k]
                             } else {
                               inik$mu <- m0_inits[k,]
                             }
                             
                             components[[k]] <<- ComponentClass$new(inik)
                           }
                           
                           a0 <<- a0
                           alpha0 <<- alpha0
                           
                           a <<- a
                           alpha <<- alpha
                           
                           log_pi <<- digamma(alpha) - digamma(sum(alpha))
                           log_rho[1] <<- digamma(a0) - digamma(2*a0)
                           log_rho[2] <<- digamma(a0) - digamma(2*a0)
                         },
                         vb_Estep_update = function() {
                            for(k in 1:K) {
                              log_q_nk[, k] <<- signals[[1]]$vb_Estep_update(X)
                              log_q_nk[, k] <<- log_q_nk + log_pi[k] + log_rho[1] 
                            }
                            log_q_nk[, K + 1] <- log_rho[2] + log(bg_density)
                             
                            # apply softmax function 
                            Z        <- apply(log_q_nk, 1, log_sum_exp)
                            log_r_nk <- log_q_nk - Z              # log of 10.49
                            r_nk     <- apply(log_r_nk, 2, exp)     # 10.49
                         },
                         vb_Mstep_update = function() {
                           # M-Step
                           for(k in 1:K) {
                             components[[k]]$vb_Mstep_update(X, r_nk[,k])
                           }
                           
                           alpha <<- alpha0 + colSums(r_nk[,1:K])
                           a[1] <<- a0 + sum(r_nk[,1:K])
                           a[2] <<-  a0 + sum(r_nk[, (K + 1)])
                         },
                         vb_EM_update = function() {
                           vb_Estep_update()
                           vb_Mstep_update()
                           
                           # update rho and we and s and z
                           update_mixing_weights()
                         },
                         compute_ELBO = function() {
                           n <- colSums(r_nk) + 1e-10 
                           
                           lb <- 0
                           for (k in 1:K) {
                             lb <- lb + signal[[1]]$components[[k]]$get_log_lb()
                           }
                           
                           lb <- lb + n[K+1] * log(bg_density) # add background
                           lb_pz  <- sum(n[1:K] %*% log_pi)    # 10.72
                           lb_pz <- lb_pz + sum(n[1:K]) * log_rho[1] + n[K+1] * log_rho[2] # add background weight
                           lb_qz  <- sum(r_nk * log_r_nk)    # 10.75
                           lb_pp  <- sum((alpha0 - 1)*log_pi) + lgamma(sum(K*alpha0)) -
                             K*sum(lgamma(alpha0))        # 10.73
                           lb_qp  <- sum((alpha - 1)*log_pi) + lgamma(sum(alpha)) -
                             sum(lgamma(alpha)) # 10.76
                           lb_pd <-lgamma(2* a0 ) - 2 * lgamma(a0) + sum(c(a0 - 1, a0 - 1) * log_rho) # p background
                           lb_qd <-lgamma(sum(a)) - sum(lgamma(a)) + sum((a-1) * log_rho)  # q background
                           
                           # Sum all parts to compute lower bound
                           L <- lb + lb_pz + lb_pp + lb_pd - lb_qz - lb_qp - lb_qd
                         },
                         update_mixing_weights = function() {
                           w <<- gtools::rdirichlet(1, alpha)
                           rho <<- rbeta(a[1], a[2])
                           
                           for(i in 1:N) {
                             s[i] <<- rbinom(1, 1, rho)
                             if(s[i] == 1) {
                               z[i] <<- sample(1:K, 1, prob = w)
                             } else {
                               z[i] <<- 0
                             }
                           }
                         },
                         get_loglik = function(x) {
                           log_probs <- rep(0, K)
                           testthat::expect_equal(sum(z != 0), sum(s == 1))
                           
                           for(k in 1:K) {
                             log_probs[k] <- log(w[k]) + components[[k]]$get_loglik(x)
                           }
                           
                           return(matrixStats::logSumExp(log_probs))
                         },
                         generate_sample = function(n){
                           cluster_allocations <- base::sample(1:K, n, replace=T, prob = w)
                           out <- matrix(0, n, D)
                           
                           # for all existing clusters
                           mu_list <- list()
                           Sigma_list <- list()
                           for(k in 1:K){
                             mu_list[[k]] <- components[[k]]$mu
                             Sigma_list[[k]] <- components[[k]]$Sigma
                             
                             subset <- (cluster_allocations == k)
                             Nk <- sum(subset) 
                             if(Nk > 0){
                               out[subset, ] <- components[[k]]$generate_data(Nk)
                             }
                           }
                           list(X = out, z = cluster_allocations, mu = mu_list, Sigma = Sigma_list)
                         }
                       ))
