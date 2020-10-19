suppressPackageStartupMessages(require(matrixcalc))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(Matrix))

# Mixture density using predictive t-distribution
mixture_pdf_t <- function(model, data){
  mixture <- vector(mode = "numeric", length = NROW(data))
  for (k in 1:length(model$nu)) {
    Lk <- solve((((model$nu[k] + 1 - NROW(model$m)) * model$kappa[k]) /
                    (1 + model$kappa[k])) * model$W[,,k])
    mixture <- mixture + (model$alpha[k]/sum(model$alpha)) *
      dmvt(x = cbind(data$x,data$y), delta = model$m[, k],
           sigma = Lk, df = model$nu[k] + 1 - NROW(model$m),
           log = FALSE, type = "shifted")
  }
  return(mixture)
}
# Use the log sum exp trick for having numeric stability
log_sum_exp <- function(x) {
  # Computes log(sum(exp(x))
  offset <- max(x)
  s <- log(sum(exp(x - offset))) + offset
  i <- which(!is.finite(s))
  if (length(i) > 0) { s[i] <- offset }
  return(s)
}

softmax <- function(...){
  logx <- unlist(list(...))
  normx <- logx - max(logx)
  exp(normx) / sum(exp(normx))
}

vb_gmm <- function(X, bg_lower =  apply(X,2, min), bg_upper =  apply(X,2, max), 
                   K = 3,
                   alpha0 = 1/K, 
                   a0 = 0.5,
                   m0 = c(colMeans(X)),
                   kappa0 = 1, 
                   nu0 = NCOL(X) + 1,
                   W0 = diag(1, NCOL(X)), 
                   m_init = NULL,
                   niters = 50,
                   epsilon_conv = 1e-6,
                   animated = FALSE) {
  # Compute logB function
  logB <- function(W, nu){
    D <- NCOL(W)
    return(-0.5*nu*log(det(W)) - (0.5*nu*D*log(2) + 0.25*D*(D - 1)*log(pi) + 
                                    sum(lgamma(0.5 * (nu + 1 - 1:D)))))  #log of B.79
  }
  
  X <- as.matrix(X)
  D <- NCOL(X)              # Number of features
  N <- NROW(X)              # Number of observations
  
  W0_inv <- solve(W0)     # Compute W^{-1}
  L <- rep(-Inf, niters)  # Store the lower bounds
  
  r_nk = log_q_nk = log_r_nk <- matrix(0, nrow = N, ncol = K + 1)
  x_bar <- matrix(0, nrow = D, ncol = K)
  S <- W <- array(0, c(D, D, K ))
  log_pi <- log_Lambda <- rep(0, K)
  log_rho <- rep(0, 2) # signal - background ratio
  
  dt_all <- data.table(x = numeric(), y = numeric(), z = numeric(),
                       iter = numeric())

  m  <- t(kmeans(X, K, nstart = 25)$centers)  # Mean of Gaussianif(!is.null(m_init)) {
  if(!is.null(m_init)) {
    for(k in 1:K) {
      if(ncol(m_init) >= k) {
        m[, k] = m_init[, k]
      }
    }
  }

  kappa <- rep(kappa0, K) # Scale of precision matrix
  nu <- rep(nu0, K) # Degrees of freedom
  alpha <- rep(alpha0, K) # Dirichlet parameter
  a <- c(a0, a0) # Beta parameter for prior ratio of signal

  log_pi  <- digamma(alpha) - digamma(sum(alpha))
  log_rho <- rep(digamma(a0) - digamma(2*a0), 2)

  for (k in 1:K) {
    W[,,k] <-  W0  # Scale matrix for Wishart
    log_Lambda[k] <- sum(digamma((nu[k] + 1 - c(1:D))/2)) +
      D*log(2) + log(det(W[,,k]))
  }

  # Iterate to find optimal parameters
  for (i in 1:niters) {
    if(i %% 5 == 0 && animated) {
      print(plot_result(X, apply(r_nk, 1, function(p) sample(1:length(p), 1, prob = p))))
    }
    
    ##-------------------------------
    # Variational E-Step
    ##-------------------------------
    for (k in 1:K) {
      Xbar <- sweep(X, 2, m[, k], "-")
      log_q_nk[, k] <- log_pi[k] + log_rho[1] - D/2 *log(2*pi) + 0.5*log_Lambda[k] -
        0.5*(D/kappa[k]) - 0.5*nu[k] * diag(Xbar %*% W[,,k] %*% t(Xbar)) # log of 10.67
    }

    log_q_nk[, K + 1] <- log_rho[2] - log(prod(bg_upper - bg_lower))
    
    # Apply softmax function & compute log_r_nk by using log_sum_exp
    Z        <- apply(log_q_nk, 1, log_sum_exp)
    log_r_nk <- log_q_nk - Z              # log of 10.49
    r_nk     <- apply(log_r_nk, 2, exp) # 10.49

    ##-------------------------------
    # Variational M-Step
    ##-------------------------------
    N <- colSums(r_nk) + 1e-10  # 10.51
    for (k in 1:K) {
      x_bar[, k] <- (r_nk[ ,k] %*% X) / N[k]   # 10.52
      x_cen <- sweep(X, 2, x_bar[, k], "-")
      S[, , k]  <- t(x_cen) %*% (x_cen * r_nk[, k]) / N[k]  # 10.53
    }
    
    # Update Dirichlet parameter
    alpha <- alpha0 + N[1:K]  # 10.58
    # # Compute expected value of mixing proportions
    pis <- (alpha0 + N[1:K]) / (K * alpha0 + sum(N[1:K]))

    # Update a parameter
    a[1] = a0 + sum(N[1:K])
    a[2] = a0 + N[K + 1]
    
    # compute expected value of delta ratio of singal to background
    rho = a[1] / (a[1] + a[2])

    # Update parameters for Gaussia-nWishart distribution
    kappa <- kappa0 + N[1:K]    # 10.60
    nu   <- nu0 + N[1:K]  # 10.63
    for (k in 1:K) {
      # 10.61
      m[, k]   <- (1/kappa[k]) * (kappa0*m0 + N[k]*x_bar[, k])
      
      # 10.62
      W[, , k] <- W0_inv + N[k] * S[,,k] +
        ((kappa0*N[k])/(kappa0 + N[k])) *
        tcrossprod((x_bar[, k] - m0))
      W[, , k] <- solve(W[, , k])
    }
    
    # Update expectations over \pi and \Lambda and \deta
    # 10.66
    log_pi <- digamma(alpha) - digamma(sum(alpha))
    log_rho = digamma(a) - digamma(sum(a))
    
    for (k in 1:K) { # 10.65
      log_Lambda[k] <- sum(digamma((nu[k] + 1 - 1:D)/2)) +
        D*log(2) + log(det(W[,,k]))
    }
    
    # Need to look at this !!!
    ##-------------------------------
    # Variational lower bound
    ##-------------------------------
    lb_px = lb_pml = lb_pml2 = lb_qml <- 0
    for (k in 1:K) {
      # 10.71
      lb_px <- lb_px + N[k] * (log_Lambda[k] - D/kappa[k] - nu[k] *
                                   matrix.trace(S[,,k] %*% W[,,k]) -
                                   nu[k]*t(x_bar[,k] -
                                   m[,k]) %*% W[,,k] %*% (x_bar[,k] - m[,k]) -
                                   D*log(2*pi) )
      # 10.74
      lb_pml <- lb_pml + D*log(kappa0/(2*pi)) + log_Lambda[k] -
        (D*kappa0)/kappa[k] - kappa0*nu[k]*t(m[,k] - m0) %*%
        W[,,k] %*% (m[,k] - m0)
      # 10.74
      lb_pml2 <- lb_pml2 + nu[k] * matrix.trace(W0_inv %*% W[,,k])
      # 10.77
      lb_qml <- lb_qml + 0.5*log_Lambda[k] + 0.5*D*log(kappa[k]/(2*pi)) -
        0.5*D - (-logB(W = W[,,k], nu = nu[k]) -
                   0.5*(nu[k] - D - 1)*log_Lambda[k] + 0.5*nu[k]*D)
    }
    lb_px  <- 0.5 * lb_px  - N[K + 1] * log(prod(bg_upper - bg_lower)) # 10.71  
    lb_pml <- 0.5*lb_pml + K*logB(W = W0,nu = nu0) + 0.5*(nu0 - D - 1) *
      sum(log_Lambda) - 0.5*lb_pml2 # 10.74
    lb_pz  <- sum(r_nk[,1:K] %*% log_pi)    # 10.72
    lb_pz <- lb_pz + sum(N[1:K]) * log_rho[1] + N[K+1] * log_rho[2] # add background weight
    lb_qz  <- sum(r_nk * log_r_nk)    # 10.75
    lb_pp  <- sum((alpha0 - 1)*log_pi) + lgamma(sum(K*alpha0)) -
      K*sum(lgamma(alpha0))        # 10.73
    lb_qp  <- sum((alpha - 1)*log_pi) + lgamma(sum(alpha)) -
      sum(lgamma(alpha)) # 10.76
    lb_prho <-lgamma(2* a0 ) - 2 * lgamma(a0) + sum(c(a0 - 1, a0 - 1) * log_rho) # p background
    lb_qrho <-lgamma(sum(a)) - sum(lgamma(a)) + sum((a-1) * log_rho)  # q background
    # Sum all parts to compute lower bound
    L[i] <- lb_px + lb_pz + lb_pp + lb_pml + lb_prho - lb_qz - lb_qp - lb_qml  - lb_qrho

    # Check if lower bound decreases
    if (i > 1 && L[i] < L[i - 1]) { message("Warning: Lower bound decreases!\n"); }
    # Check for convergence
    if (i > 1 && abs(L[i] - L[i - 1]) < epsilon_conv) { break }
    # Check if VB converged in the given maximum iterations
    if (i == niters) {warning("VB did not converge!\n")}
  }


  # Rearranging the order so that cluster 1 has smallest mu_x
  sort_order <- order(m[1,])
  pis <- pis[sort_order]
  m <- m[, sort_order]
  W <- W[,, sort_order]
  kappa <- kappa[sort_order]
  nu <- nu[sort_order]
  alpha <- alpha[sort_order]
  r_nk <- r_nk[, c(sort_order, K +1)]

  obj <- structure(list(X = X, K = K, N = N, D = D, pi = pis, rho = rho,
                        a = a, alpha = alpha, r_nk = r_nk,  m = m, W = W,
                        kappa = kappa, nu = nu, L = L[2:i], bg_lower = bg_lower,
                        bg_upper = bg_upper, dt_all = dt_all), class = "vb_gmm")
  return(obj)
}