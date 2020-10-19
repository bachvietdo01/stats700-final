EMFMMSignal <- setRefClass("EMFMMSignal",
                       fields = list(
                         components = "list",
                         K = "integer",
                         D = "integer",
                         X = "matrix",
                         mz_hat = "matrix",
                         ComponentClass ="refObjectGenerator",
                         init_pars = "list",
                         N = "integer",
                         bg_density = "numeric",
                         rho = "numeric",
                         w = "numeric",
                         s = "numeric",
                         z = "numeric"),

                       methods = list(
                         initialize = function(K, D, X, s, z, 
                                               rho = 0.5, w = rep(1/K, K),
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
                           
                           rho <<- rho
                           w <<- w
                           bg_density <<- bg_density
                           mz_hat <<- matrix(1e-10, nrow = N, ncol = (K + 1)) # N x (K + 1) where 1st col is z_i == 0
                           
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
                         },
                         em_update = function() {
                           # E-Step
                           for(i in 1:N) {
                             x <- X[i,]
                             mz_hat[i, 1] <<- log(1 - rho) + log(bg_density)
                             
                             for(k in 1: K) {
                               mz_hat[i, k + 1] <<- log(rho) + log(w[k]) + components[[k]]$get_loglik(x)
                             }
                             
                             mz_hat[i,] <<- softmax(mz_hat[i,])
                             if(any(is.na(mz_hat))) {
                               stop("update mz_hat")
                             }
                             z[i] <<- base::sample(1 : (K+1), 1, prob = mz_hat[i,]) - 1
                           }
                           
                           # update s
                           s[z == 0] <<- 0
                           s[z != 0] <<- 1
                           
                           # M-Step
                           for(k in 1:K) {
                             # 1st column zi == 0
                             # update mu, Sigma
                             zk_hat <- mz_hat[, k + 1]
                             if(sum(zk_hat) > 0) {
                               components[[k]]$M_step_update(X, zk_hat)
                             }
                             
                             # update w
                             w[k] <<-  log(sum(mz_hat[, k + 1]))
                           }
                           
                           w <<- softmax(w)
                           rho <<- sum(mz_hat[,2:(K+1)]) / N
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
