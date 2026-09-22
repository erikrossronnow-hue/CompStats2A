# can write the gaussian stuff as a function here, test the function in the sandbox
# with a new script please. 


# mu
mu <- uniroot(
  sum(x * z) - sum(x * exp(y * x)),
  interval = c(0, 10),
  x = x,
  z = z
)$root
  
# sigma
sigma <- 1 / sqrt(sum(x^2 * exp(mu * x)))

# f function
log_f <- function(y, x, z) {
  sum0 <- sum(y*x*z)
  sum1 <- sum(exp(y*x))
  sum0 - sum1
} # Need to add protection for y <0 ???

# proposal function 
log_q <- function(y, mu, sigma) {
  dnorm(y, mean = mu, sd = sigma, log = TRUE)
}


# Find M
# Not quite sure how but probably also need to be some kind of maximizing function

# Recall that the shape parameter has to fulfill: r >= 1
gaussian_random <- function(N, mu, sigma, M, x, z) {
  
  y <- rnorm(N, mean = mu, sd = sigma)
  u <- runif(N)
  
  accept <- y >= 0 & u <= exp(log_f(y, x, z) -
                                     log_q(y, mu, sigma) -
                                     log(M))
  
  y[accept]
}
