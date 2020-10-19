gen_mix_gmm_and_unif <- function(N, Umin, Umax, mu, Sigma, rho, pi) {
  d = length(Umin)
  
  X = matrix(rep(0, d* N), ncol = d)
  s = rep(0, N)
  z = rep(0, N)
  
  for (n in 1: N) {
    sn = rbinom(1, 1, rho) # sample s
    s[n] = sn
    
    if (sn == 0) { # sample background from a uniform distribution
      for (p in 1 :d) {
        X[n,p] = runif(1, min = Umin[p], max = Umax[p])
      } 
    } else { # sample signal from a mixture of k multivariate normal
      zn = sample(1:K, 1, prob = pi)
      z[n] = zn
      # Get mu and Sigma for kth normal component
      muk = mu[, zn]
      Sigmak = Sigma[, , zn]
      x = rep(0, d)  
      while(TRUE) {
        x = mvtnorm::rmvnorm(1, mean = muk, sigma = Sigmak)
        for(i in 1 :d) {
          if (x[i] < Umin[i] || x[i] > Umax) { #outside of allowed region
              next
          }
        }
        
        X[n,] = x
        break
      }
    }
  }
  
  return(list(X = X, s = s, z = z))
}  


# Construct initial parameters
N = 1000;
d = 2;
K = 2;

# Construct mu and Sigma
I = diag(d);
mu0 = c(-10, 10);

mu = matrix(rep(0, d * K), nrow = 2)
mu[,1] = c(-5, 0)
mu[,2] = c(1, 0)

sigma.2 = 3; # know varaince
v = c();
for (k in 1: K) {
  v = c(v, as.vector(sigma.2 * I))
}
Sigma = array(v, c(d, d, K))

# Construct Umin and Umax
Umin = c()
Umax = c()
for (p in 1:d) {
  Umax = c(Umax, max(mu[p,]) + 3 * sigma.2) 
  Umin = c(Umin, min(mu[p,]) - 3 * sigma.2);
}

# Set rho and pi
rho = 0.75;
pi = c(0.7, 0.3)

data = list()
output = gen_mix_gmm_and_unif(N, Umin, Umax, mu, Sigma, rho, pi)
data$X = output$X
data$s = output$s
data$z = output$z
data$mu = mu
data$rho = rho
data$pi = pi
data$Sigma = Sigma
data$lower_bg = Umin
data$upper_bg = Umax
data$N = N
data$K = K
data$d = d
rm(list=setdiff(ls(), "data"))

df = data.frame(x = data$X[, 1], y = data$X[, 2], z = data$z)
ggplot(df, aes(x, y, col = z)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")