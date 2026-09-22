# can write the gaussian stuff as a function here, test the function in the sandbox
# with a new script please. 

# For the standard gaussian 

log_f <- function(y, x, z) {
  sum0 <- sum(y*x*z)
  sum1 <- sum(exp(y*x))
  sum0 - sum1
}

log_q <- function(y, mu, sigma) {
  dnorm(y, mean = mu, sd = sigma, log = TRUE)
}

gaussian_random <- function(N, mu, sigma, M, x, z) {
  
  y <- rnorm(N)
  u <- runif(N)
  
  accept <- y >= 0 & u <= exp(log_f(y, x, z) -
                                log_q(y, 0, 1))
  y[accept]
}



