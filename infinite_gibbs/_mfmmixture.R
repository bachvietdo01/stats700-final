MFMMixture <- setRefClass("MFMMixture",
                       fields = list(
                         sig_components = "list",
                         K = "integer",
                         D = "integer",
                         X = "matrix",
                         init_pars = "list",
                         N_total = "integer",
                         s = "numeric",
                         z = "numeric",
                         rho = "numeric",
                         w = "numeric",
                         kappa0 = "numeric",
                         nu0 = "numeric",
                         m0 = "numeric",
                         L0 = "matrix",
                         S0 = "matrix",
                         m0_init = "matrix",
                         alpha = "numeric",
                         alpha0 = "numeric",
                         lower_bg = "numeric",
                         upper_bg = "numeric",
                         log_V_n = "numeric"),

                       methods = list(
                         initialize = function(K, D, X, z, s,
                                               kappa0 = 1.0, 
                                               nu0 = D+2, 
                                               m0 = rep(0, D), 
                                               S0 = diag(rep(1, D)), 
                                               alpha = 0.1,
                                               alpha0 = 0.5,
                                               m0_init = matrix(0, ncol = K,
                                                                nrow = D)){
                           K <<- as.integer(K)
                           D <<- as.integer(D)
                           N_total <<- nrow(X)
                           s <<- s
                           z <<- z
                           X <<- X
                           alpha0 <<- alpha0
                           alpha <<- alpha

                           kappa0 <<- kappa0
                           nu0 <<- nu0
                           m0 <<- m0
                           S0 <<- S0
                           L0 <<- chol(S0 + kappa0 * outer(m0, m0))
                           
                           for(k in 1:K){
                             if(ncol(m0_init) >= k) {
                               m0k <- m0_init[,k]
                             } else {
                               m0k <- m0
                             }
                             
                             L0k <- chol(S0 + kappa0 * outer(m0k, m0k))
                             init_pars <<- list(D = D, kappa0 = kappa0, nu0 = nu0, m0 = m0k, S0 = S0, L0 = L0k)
                             sig_components[[k]] <<- SignalComponent$new(init_pars, X = X[z == k, , drop=FALSE])
                            }
                           
                           rho <<- rbeta(1, alpha0, alpha0)
                           w <<- rep(1/K, K)
                           lower_bg <<- apply(X, 2, min)
                           upper_bg <<- apply(X, 2, max)
                           
                           log_V_n <<- calculate_log_V_n(alpha, sum(s == 1), 100)
                         },
                         update_X = function(X){
                           # given a new X value, but conditioning on the current clustering z,
                           # we need to update m and L for each component
                           # (the rest will remain constant)
                           X <<- X
                           for(k in 1:K){
                             X_k <- X[z == k, , drop=FALSE]
                             N_k <- sum(z == k)
                             kappa_k <- (kappa0 + N_k)
                             m_k <- (kappa0 * m0 + colSums(X_k)) / kappa_k
                             S_k <- S0 + kappa0 * outer(m0, m0) + t(X_k) %*% X_k
                             sig_components[[k]]$m <<- m_k
                             sig_components[[k]]$L <<- chol(S_k - kappa_k * outer(m_k, m_k))
                             sig_components[[k]]$S <<- S_k
                           }
                         },
                         update_rho = function() {
                           N1 <- sum(s == 1)
                           rho <<- rbeta(1, alpha0 + N1, alpha0 + N_total - N1)
                         },
                         update_w = function() {
                           cnt_z <- rep(0, K)
                           for(k in 1:K) {
                             cnt_z[k] <- sum(z == k)
                           }
                           w <<- as.vector(MCMCpack::rdirichlet(1, rep(alpha, K) + cnt_z))
                           
                         },
                         add_background = function(i) {
                           s_probs <- rep(0, 2)
                           s_probs[2] <- log(1 - rho) - log(prod(upper_bg - lower_bg))
                           
                           log_sig_probs <- rep(0, K)
                           N_s <- sum(z != 0)
                           testthat::expect_equal(N_s, sum(s == 1))
                           
                           for(k in 1:K) {
                             w_k <- w[k]
                             log_sig_probs[k] <- log(rho) + log(w_k) +  
                               sig_components[[k]]$get_loglik(X[i, ])
                           }
                           
                           s_probs[1] <- matrixStats::logSumExp(log_sig_probs)
                           
                           if(! is.na(s_probs[1])) {
                             s_probs = softmax(s_probs)
                           } else {
                             s_probs= c(0,1)
                           }
                           
                           s[i] <<- rbinom(1, 1, s_probs[1])
                           if (s[i] == 0) {
                             z[i] <<- 0
                           }
                         },
                         rm_background = function(i) {
                           s[i] <<- 0
                           if(z[i] != 0) {
                             rm_signal(i)
                             z[i] <<- 0
                           }
                         },
                         add_signal = function(i, k){
                           x <- X[i, ]
                           z[i] <<- k
                           if(k > K){
                             new_sig_component()
                           }
                           
                           sig_components[[k]]$add_sample(x)
                         },
                         rm_signal = function(i){
                           x <- X[i, ]
                           k <- z[i]
                           
                           sig_components[[k]]$rm_sample(x)
                           
                           if(sig_components[[k]]$is_empty()){
                             rm_sig_component(k)
                           }
                         },
                         new_sig_component = function(which_ind = NULL){
                           K <<- K + 1L
                           if(is.null(which_ind)){
                             sig_components[[K]] <<- SignalComponent$new(init_pars)
                           }else{
                             sig_components[[K]] <<- SignalComponent$new(init_pars, X = X[which_ind, , drop=FALSE])
                             z[which_ind] <<- K
                           }
                         },
                         rm_sig_component = function(k){
                           sig_components[[k]] <<- NULL
                           K <<- K - 1L
                           # update the cluster numbers
                           z <<- ifelse(z > k, z-1L, z)
                         },
                         signal_collapsed_gibbs_for_obs_i = function(i){
                           x <- X[i, ]
                           logprobs <- rep(0, K+1)
                           # for existing clusters
                           for(k in 1:K){
                             logprior <- log(sum(z[-i] == k) + alpha)
                             loglik <- sig_components[[k]]$posterior_predictive(x)
                             logprobs[k] <- logprior + loglik
                             # cat("k", k, "logprior", logprior, "loglik", loglik, "\n")
                           }
                           # cat("logprobs" ,logprobs, "\n")
                           # for a new cluster
                           Ltemp <- add_sample_and_get_cholS_helper(x, m0, kappa0, nu0, L0)
                           
                           # temp <- helper_likelihood(1, D, kappa0+1, nu0+1, Ltemp) - helper_likelihood(0, D, kappa0, nu0, L0)
                           temp <- -0.5*D*log(pi) - 0.5*D*log((kappa0+1)/kappa0) - 0.5*(nu0+1)*logdet_chol(Ltemp) + lgamma(0.5*(nu0 + 1)) +
                             0.5*nu0*logdet_chol(L0) - lgamma(0.5*(nu0 + 1 - D))
                           logprobs[K+1] <- log(alpha) + log_V_n[K+1] - log_V_n[K] + temp
                           
                           k0 = base::sample(1:(K+1), 1, prob = softmax(logprobs))
                           k0
                         },
                         collapsed_gibbs = function(){
                           for(i in 1:N_total){
                             # add backgroumd
                             rm_background(i)
                             add_background(i)
                             update_rho()
                             
                             log_V_n <<- calculate_log_V_n(alpha, sum(s == 1), 100)
                             
                             # add signal
                             if(s[i] == 1) {
                               k0 <- signal_collapsed_gibbs_for_obs_i(i)
                               add_signal(i, k0)
                               update_w()
                             }
                             
                             testthat::expect_equal(sum(z != 0), sum(s == 1))
                           }
                         },
                         merge_split = function(){
                           points <- 1:N_total
                           signals <- points[z != 0]
                           
                           if (length(signals) > 0) {
                             i <- sample(signals, 1)
                             j <- sample(setdiff(signals, i), 1)
                             testthat::expect_true(z[i] != 0)
                             testthat::expect_true(z[j] != 0)
                           
                             if(z[i] == z[j]){
                               propose_split(i, j, sig_components[[z[i]]])
                               testthat::expect_equal(sum(s == 1), sum(z != 0))
                             } else{
                               propose_merge(i, j)
                               testthat::expect_equal(sum(s == 1), sum(z != 0))
                             }
                             update_w()
                           }
                         },
                         propose_split = function(i, j, S_current){
                           S_i <- SignalComponent$new(init_pars)
                           S_j <- SignalComponent$new(init_pars)
                           S_i$add_sample(X[i, ])
                           S_j$add_sample(X[j, ])
                           S_ind <- which(z %in% c(z[i], z[j]))
                           # temp_z <- rep(0L, length(S_ind))
                           # temp_z[S_ind == i] <- 1L
                           # temp_z[S_ind == j] <- 2L
                           temp_z <- rep(0L, N_total)
                           temp_z[c(i, j)] <- c(1L, 2L)
                           MH_logratio <- 0
                           if(length(S_ind) > 2){
                             ind_perm <- sample((S_ind))
                             # temp cluster allocations within S_ind
                             for(kk in setdiff(ind_perm, c(i, j))){
                               # kk <- S_ind[k]
                               x <- X[kk, ]
                               # choose whether add observation k to S_i or S_j
                               p_i <- S_i$N * exp(S_i$posterior_predictive(x))
                               p_j <- S_j$N * exp(S_j$posterior_predictive(x))
                               prob_i <- p_i / (p_i + p_j)
                               if(runif(1) < prob_i){
                                 S_i$add_sample(x)
                                 temp_z[kk] <- 1L
                                 MH_logratio <- MH_logratio + log(prob_i)
                               } else{
                                 S_j$add_sample(x)
                                 temp_z[kk] <- 2L
                                 MH_logratio <- MH_logratio + log(1-prob_i)
                               }
                             }
                           }
                           logprob_proposed <- S_j$marginal_loglik() + S_i$marginal_loglik()
                           logprob_current <- S_current$marginal_loglik()
                           MH_logratio <- - MH_logratio + logprob_proposed - logprob_current + log(alpha) + lgamma(S_i$N) + lgamma(S_j$N) - lgamma(S_i$N + S_j$N)
                           # accept or reject the constructed proposal
                           if(runif(1) < exp(MH_logratio)){
                             # cat("accepted split cluster", z[i], "with prob", exp(MH_logratio), "\n")
                             rm_sig_component(z[i])
                             # cat("splitting elements", S_ind, "into", S_ind[temp_z == 1L], "and", S_ind[temp_z == 2L], "\n")
                             # new_component(which_ind = which(temp_z == 1L))
                             sig_components[[K+1]] <<- S_i
                             sig_components[[K+2]] <<- S_j
                             z[temp_z == 1L] <<- K+1L
                             z[temp_z == 2L] <<- K+2L
                             K <<- K + 2L
                             # new_component(which_ind = which(temp_z == 2L))
                             # cat("new splits:", which(temp_z == 1L), "and", which(temp_z == 2L), "\n")
                           }
                         },
                         propose_merge = function(i, j){
                           S_ind <- which(z %in% c(z[i], z[j]))
                           S_merged <- SignalComponent$new(init_pars, X = X[S_ind, , drop=FALSE])
                           S_i <- SignalComponent$new(init_pars)
                           S_j <- SignalComponent$new(init_pars)
                           S_i$add_sample(X[i, ])
                           S_j$add_sample(X[j, ])
                           MH_logratio <- 0
                           if(length(S_ind) > 2){
                             # imaginary clusters S_i and S_j
                             ind_perm <- sample(setdiff(S_ind, c(i, j)))
                             for(k in 1:length(ind_perm)){
                               kk <- ind_perm[k]
                               x <- X[kk, ]
                               # choose whether add observation k to S_i or S_j
                               p_i <- S_i$N * exp(S_i$posterior_predictive(x))
                               p_j <- S_j$N * exp(S_j$posterior_predictive(x))
                               prob_i <- p_i / (p_i + p_j)
                               if(z[kk] == z[i]){
                                 S_i$add_sample(x)
                                 MH_logratio <- MH_logratio + log(prob_i)
                               } else{
                                 S_j$add_sample(x)
                                 MH_logratio <- MH_logratio + log(1-prob_i)
                               }
                             }
                           }
                           logprob_current <- S_j$marginal_loglik() + S_i$marginal_loglik()
                           logprob_proposed <- S_merged$marginal_loglik()
                           MH_logratio <- MH_logratio + logprob_proposed - logprob_current - log(alpha) - lgamma(S_i$N) - lgamma(S_j$N) + lgamma(S_i$N + S_j$N)
                           if(runif(1) < exp(MH_logratio)){
                             # cat("accepted merge with prob", exp(MH_logratio), "\n")
                             rm_sig_component(z[i])
                             rm_sig_component(z[j])
                             new_sig_component(which_ind = S_ind)
                           }
                           rm(S_i, S_j, S_merged)
                         },
                         generate_sample = function(n){
                           cluster_probs <- sapply(sig_components, function(x)x$N)
                           cluster_allocations <- base::sample(1:K, n, replace=T, prob = cluster_probs)
                           out <- matrix(0, n, D)
                           # for all existing clusters
                           mu_list <- list()
                           Sigma_list <- list()
                           for(k in 1:K){
                             Sigma_list[[k]] <- sig_components[[k]]$Sigma
                             mu_list[[k]] <- sig_components[[k]]$mu
                             subset <- (cluster_allocations == k)
                             if(sum(subset) > 0){
                               out[subset, ] <- mvtnorm::rmvnorm(sum(subset), mean = mu_list[[k]], sigma = Sigma_list[[k]])
                             }
                           }

                           N1 <- sum(s == 1)
                           rhos <- rbeta(1, alpha0 + N1, alpha0 + N_total - N1)
                           # # for the new cluster
                           # subset <- (cluster_allocations == k+1)
                           # if(sum(subset) > 0){
                           #   Sigma <- solve(rWishart(1, nu0, solve(S0))[, , 1])
                           #   mu <- as.numeric(mvtnorm::rmvnorm(1, mean = m0, sigma = 1/kappa0 * Sigma))
                           #   cat("Started a new cluster, with mean", mu, "and covariance ", Sigma, "\n")
                           #   out[subset, ] <- mvtnorm::rmvnorm(sum(subset), mean = mu, sigma = Sigma)
                           # }
                           list(X = out, z = cluster_allocations, mu = mu_list, Sigma = Sigma_list, rho = rhos)
                         },
                         tidy_up_clusters = function() {
                           testthat::expect_equal(sum(z != 0), sum(s == 1))
                           zs <- z[z != 0]
                           cstar <- unique(zs)
                           
                           # store index of old cluster numbers
                           old_sig_comps <- sig_components
                           idx = list()
                           for(c in cstar) {
                             idx[[c]] <- (z == c)
                           }
                           
                           # reassign cluster numbers
                           sig_components <<- list()
                           K <<- length(cstar)
                           k <- 1
                           for(c in cstar) {
                             sig_components[[k]] <<- old_sig_comps[[c]]
                             z[idx[[c]]] <<- rep(k, sum(idx[[c]]))
                             k = k + 1
                           }
                           testthat::expect_equal(sum(z != 0), sum(s == 1))
                           update_w()
                         },
                         check_L = function(){
                           for(k in 1:K){
                             X_k <- X[z == k, , drop=FALSE]
                             N_k <- nrow(X_k)
                             kappa_k <- kappa0 + N_k
                             m_k <- (kappa0*m0 + colSums(X_k)) / kappa_k
                             S_k <- S0 + t(X_k) %*% X_k
                             L_k <- chol(S_k - kappa_k * outer(m_k, m_k))
                             testthat::expect_equal(L_k, sig_components[[k]]$L)
                           }
                         },
                         plot = function(){
                           # r1 <- extendrange(X[, 1])
                           # r2 <- extendrange(X[, 2])
                           # xval <- seq(r1[1], r1[2], 0.02)
                           # yval <- seq(r2[1], r2[2], 0.02)
                           # dftemp <- melt(outer(xval , yval))
                           # x <- xval[df.temp$Var1]
                           # y <- yval[df.temp$Var2]
                           # mat <- cbind(x, y)
                           # df <- data.frame(c())
                           mu_list <- list()
                           Sigma_list <- list()
                           
                           for(k in 1:K){
                             Sigma_list[[k]] <- sig_components[[k]]$Sigma
                             mu_list[[k]] <- sig_components[[k]]$mu
                           }
                           plot_mixture(mu_list, Sigma_list, X, z)
                         }
                       ))
