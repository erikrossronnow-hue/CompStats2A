poisson <- read.csv(
  "https://math-ku.github.io/compstat/data/poisson.csv"
)

z <- poisson$z
x <- poisson$x

log_f <- function(y) {
  sum(y * z * x - exp(y * x))
}


mu <- 0.24
sigma <- 4e-2
log_g <- function(y) {
  dnorm(y, mean = mu, sd = sigma, log = TRUE)
}

y_grid <- seq(0, 1, length.out = 1000)

# log_f evaluates one y at a time
f_values <- vapply(y_grid, log_f, numeric(1))
g_values <- log_g(y_grid)

quartz(bg = "white", width = 8, height = 5)

plot(
  y_grid, f_values,
  type = "l",
  col = "blue",
  lwd = 2,
  xlim = c(0, 1),
  ylim = range(c(f_values, g_values)),
  xlab = "y",
  ylab = "Log density"
)

lines(
  y_grid, g_values,
  col = "red",
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c("log_f(y)", "log_g(y)"),
  col = c("blue", "red"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

log_ratio <- function(y) {
  log_g(y) - log_f(y)
}

result <- optimize(
  f = log_ratio, interval = c(0,200)
)

result
log_alpha <- result$objective
log_alpha
get_one_sample <- function(log_alpha){
  rejected <- 0
  repeat{
    Y <- rnorm(1, mean = mu, sd = sigma)
    U <- runif(1)

    if(Y>=0 && log(U)<=log_alpha + log_f(Y) - log_g(Y)){
      break
    }
    rejected <- rejected + 1
  }

  list(sample = Y, rejected = rejected)
}

rejection_sample <- function(N, log_alpha){
  results <- replicate(N, get_one_sample(log_alpha), simplify = FALSE)
  list(
    samples = sapply(results, function(a) a$sample),
    rejected = sum(sapply(results, function(a) a$rejected))
  )
}


result_samples <- rejection_sample(
  N = 3000,
  log_alpha = log_alpha
)

samples <- result_samples$samples
result_samples$rejected
3000 / (3000 + result_samples$rejected)

hist(
  samples,
  breaks = 20,
  probability = TRUE,
  col = "lightblue",
  border = "white",
  main = "Rejection Sampling",
  xlab = "y",
  plot = TRUE
)
