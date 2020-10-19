Mixture <- setRefClass("Mixture",
                       fields = list(
                         background = "list",
                         signal = "list",
                         Ks = "integer",
                         Kb = "integer",
                         D = "integer",
                         X = "matrix",
                         N_total = "integer",
                         s = "numeric",
                         sz = "numeric",
                         bz = "numeric",
                         rho = "numeric",
                         alpha0 = "numeric",
                         sig_alpha = "numeric",
                         bg_alpha = "numeric",
                         sig_componentClass ="refObjectGenerator",
                         sig_init_pars = "list",
                         bg_init_pars = "list",
                         lower_bg = "numeric",
                         upper_bg = "numeric",
                         type = 'character'),

                       methods = list(
                         initialize = function(Ks, Kb, D, X, sz, bz, s,
                                               sig_componentClass,
                                               sig_init_pars,
                                               bg_init_pars = NULL,
                                               alpha0 = 0.5,
                                               sig_alpha = 1 / Ks,
                                               type = 'FMM'){
                           
                           type <<- type
                           Ks <<- as.integer(Ks)
                           Kb <<- as.integer(Kb)
                           D <<- as.integer(D)
                           N_total <<- nrow(X)
                           s <<- s
                           bz <<- bz
                           sz <<- sz
                           
                           X <<- X
                           
                           alpha0 <<- alpha0
                           sig_alpha <<- sig_alpha
                           bg_alpha <<- bg_alpha
                           
                           sig_componentClass <<- sig_componentClass
                           
                           if(type == 'DPSignal') {
                             signal[[1]] <<- DPSignal$new(K = Ks, D = D, X = X, s = s, 
                                                      z = sz, alpha = sig_alpha, sig_init_pars)
                           } else if (type == 'MFMSignal') {
                             signal[[1]] <<- MFMSignal$new(K = Ks, D = D, X = X, s = s, 
                                                          z = sz, alpha = sig_alpha, sig_init_pars)
                           } else {
                             signal[[1]] <<- FMMSignal$new(K = Ks, X = X, s = s, z = sz,
                                                                   alpha = sig_alpha,
                                                                   ComponentClass = sig_componentClass,
                                                                   init_pars = sig_init_pars)
                           }
                           
                           # if no background models specified use uniform background
                           if(is.null(bg_init_pars)) {
                             lower_bg <<- apply(X, 2, min)
                             upper_bg <<- apply(X, 2, max)
                           }
                           
                           rho <<- rbeta(1, alpha0, alpha0)
                         },
                         update_rho = function() {
                           N1 <- sum(s == 1)
                           rho <<- rbeta(1, alpha0 + N1, alpha0 + N_total - N1)
                         },
                         gibbs = function(){
                           for(i in 1:N_total){
                             xi <- X[i,]
                             
                             s_probs <- rep(0, 2)
                             s_probs[1] <- log(rho) + signal[[1]]$get_loglik(xi)
                             s_probs[2] <- log(1 - rho) - log(prod(upper_bg - lower_bg))
                             
                             if(! is.na(s_probs[1])) {
                               s_probs = softmax(s_probs)
                             } else {
                               s_probs= c(0,1)
                             }
                             
                             s[i] <<- rbinom(1, 1, s_probs[1])
                           }
                          
                          update_rho()

                          signal[[1]]$set_s(s)
                          signal[[1]]$gibbs()
                          sz <<- signal[[1]]$z
                          
                          # split-merge steps for DP & MFM
                          if(type != 'FMM') {
                            signal[[1]]$merge_split()
                            signal[[1]]$tidy_up_clusters()
                          }
                        },
                        generate_sample = function(n) {
                          samp <- signal[[1]]$generate_sample(n)
                          samp$rho <- rho
                          return(samp)
                        }
                       ))
