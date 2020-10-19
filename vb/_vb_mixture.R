EMMixSigAndBg <- setRefClass("EMMixSigAndBg",
                       fields = list(
                         signal = "list",
                         Ks = "integer",
                         N = "integer",
                         D = "integer",
                         X = "matrix",
                         s = "numeric",
                         z = "numeric",
                         rho = "numeric",
                         sig_componentClass ="refObjectGenerator",
                         sig_init_pars = "list",
                         lower_bg = "numeric",
                         upper_bg = "numeric",
                         type = 'character'),

                       methods = list(
                         initialize = function(Ks, D, X, s, z,
                                               alpha0 = 1/ Ks,
                                               a0 = 0.5,
                                               sig_componentClass,
                                               sig_init_pars){
                           
                           Ks <<- as.integer(Ks)
                           D <<- as.integer(D)
                           N <<- nrow(X)
                           s <<- s
                           z <<- z
                           X <<- X
                           rho <<- rho
                           
                           sig_componentClass <<- sig_componentClass
                           lower_bg <<- apply(X, 2, min)
                           upper_bg <<- apply(X, 2, max)
                           bg_density <- 1.0 / prod(upper_bg - lower_bg)
                           
                           signal[[1]] <<- VBFMMSignal$new(K = Ks, D = D, X = X, s = s, 
                                                      z = z, alpha0 = alpha0, a0 = a0,
                                                      bg_density = bg_density, 
                                                      ComponentClass = sig_componentClass,
                                                      init_pars = sig_init_pars)
                         },
                         vb_update = function() {
                           signal[[1]]$vb_EM_update()
                           
                           rho <<- signal[[1]]$rho
                           s <<- signal[[1]]$s
                           z <<- signal[[1]]$z
                         },
                         get_loglik = function() {
                           loglik <- 0
                           
                           for(i in 1:N) {
                             xi <- X[i,]
                             loglik <- loglik + log(1 - rho) - log(prod(upper_bg - lower_bg)) +
                               signal[[1]]$get_loglik(xi)
                           }
                           
                           return(loglik)
                         }
                       ))
