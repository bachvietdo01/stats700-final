FMMSignal <- setRefClass("FMMSignal",
                       fields = list(
                         components = "list",
                         K = "integer",
                         X = "matrix",
                         D = "integer",
                         alpha = "numeric",
                         ComponentClass ="refObjectGenerator",
                         init_pars = "list",
                         N_total = "integer",
                         s = "numeric",
                         z = "numeric",
                         w = "numeric"),

                       methods = list(
                         initialize = function(K, X, s, z, 
                                               alpha = 1.0 /K, 
                                               ComponentClass = NULL,
                                               init_pars = NULL){
                           testthat::expect_true(!is.null(ComponentClass))
                           
                           ComponentClass <<- ComponentClass
                           
                           D <<- as.integer(ncol(X))
                           K <<- as.integer(K)
                           N_total <<- nrow(X)
                           s <<- s
                           z <<- z
                           X <<- X
                           
                           m0_inits <- kmeans(x = X, K)$centers
                           
                           for(k in 1:K){
                             inik <- init_pars
                             
                             if(!is.null(init_pars$m_init)) {
                               inik$m0 <- init_pars$m_init[,k]
                             } else {
                               inik$m0 <- m0_inits[k,]
                             }
                             
                             components[[k]] <<- ComponentClass$new(inik, X = X[z == k, , drop=FALSE])
                           }
                           
                           alpha <<- alpha
                           w <<- rep(alpha, K)
                         },
                         update_w = function() {
                           cnt_z <- rep(0, K)
                           for(k in 1:K) {
                             cnt_z[k] <- sum(z == k)
                           }
                           w <<- as.vector(MCMCpack::rdirichlet(1, rep(alpha, K) + cnt_z))
                         },
                         set_s = function(s) {
                           s <<- s
                         },
                         add_signal = function(i) {
                           s[i] <<- 1
                           k <- z[i]
                           
                           testthat::expect_true(k > 0)
                           components[[k]]$add_sample(X[i,])
                         },
                         rm_signal = function(i) {
                           k <- z[i]
                           s[i] <<- 0
                           z[i] <<- 0
                           
                           testthat::expect_true(k > 0)
                           components[[k]]$rm_sample(X[i,])
                         },
                         gibbs_for_obs_i = function(i){
                           # remove point from old cluster
                           if(z[i] > 0) {
                             rm_signal(i)
                           }
                           
                           # assign new cluster
                           x <- X[i, ]
                           logprobs <- rep(0, K)
                           # for existing clusters
                           for(k in 1:K){
                             logprobs[k] <- log(w[k]) + components[[k]]$get_loglik(x)
                           }
                           z[i]  <<- base::sample(1:K, 1, prob = softmax(logprobs))
                           # add point to new cluster
                           add_signal(i)
                         },
                         gibbs = function(){
                           for(i in 1:N_total){
                             if(s[i] == 1) {
                               gibbs_for_obs_i(i)
                             } else {
                               # remove from consideration
                               if(z[i] > 0) {
                                 rm_signal(i)
                               }
                             }
                             
                             update_w()
                           }
                           testthat::expect_equal(sum(z != 0), sum(s == 1))
                         },
                         get_loglik = function(x) {
                           log_sig_probs <- rep(0, K)
                           testthat::expect_equal(sum(z != 0), sum(s == 1))
                           
                           for(k in 1:K) {
                             log_sig_probs[k] <- log(w[k]) + components[[k]]$get_loglik(x)
                           }
                           
                           return(matrixStats::logSumExp(log_sig_probs))
                         },
                         generate_sample = function(n){
                           cluster_probs <- sapply(components, function(x)x$N)
                           cluster_allocations <- base::sample(1:K, n, replace=T, prob = cluster_probs)
                           out <- matrix(0, n, D)
                           # for all existing clusters
                           mu_list <- list()
                           Sigma_list <- list()
                           for(k in 1:K){
                             Sigma_list[[k]] <- components[[k]]$Sigma
                             mu_list[[k]] <- components[[k]]$mu
                             subset <- (cluster_allocations == k)
                             if(sum(subset) > 0){
                               out[subset, ] <- mvtnorm::rmvnorm(sum(subset), mean = mu_list[[k]], sigma = Sigma_list[[k]])
                             }
                           }
                           list(X = out, z = cluster_allocations, mu = mu_list, 
                                Sigma = Sigma_list)
                         }
                       ))
