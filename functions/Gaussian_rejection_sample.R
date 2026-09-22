# can write the gaussian stuff as a function here, test the function in the sandbox
# with a new script please. 
Gaussian_rejection_sample = function() {
  return(0)
}


# mu
mu <- # need to do something with roots
  
  # sigma
  sigma <- 1 / sqrt(sum(x^2 * exp(mu * x)))

# f function
log_f <- function(y, x, z) {
  sum0 <- sum(y*x*z)
  sum1 <- sum(exp(y*x))
  sum0 - sum1
}

# proposal function 
log_q <- function(y, mu, sigma) {
  dnorm(y, mean = mu, sd = sigma, log = TRUE)
}


# Find M


# Recall that the shape parameter has to fulfill: r >= 1
gaussian_random <- function(N, mu, sigma, M, x, z) {
  
  y <- rnorm(N, mean = mu, sd = sigma)
  u <- runif(N)
  
  log_accept <- log_f(y, x, z) -
    log_q(y, mu, sigma) -
    log(M)
  
  accept <- y >= 0 & log(u) <= log_accept
  
  y[accept]
}
