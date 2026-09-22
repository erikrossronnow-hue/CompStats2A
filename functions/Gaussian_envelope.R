poisson <- read.csv(
  "https://math-ku.github.io/compstat/data/poisson.csv"
)

z <- poisson$z
x <- poisson$x

log_f <- function(y) {
  sum(y * z * x - exp(y * x))
}

log_g <- function(y){
  dnorm(y, mean=0, sd=1, log=TRUE)
}

log_ratio <- function(y) {
  log_g(y) - log_f(y)
}

result <- optimize(
  f = log_ratio, interval = c(0,200)
)

log_alpha <- result$objective

get_one_sample <- function(log_alpha){
  rejected <- 0
  repeat{
    Y <- rnorm(1)
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
  xlab = "y"
)

sample_mean <- mean(samples)
sample_variance <- var(samples)

sample_mean
sample_variance
