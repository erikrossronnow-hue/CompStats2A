# log_concave simulation function using rejection sampling from chapter 6 in the book. 
# erik: we could further optimize this function by not using matrices and more vectorization
# or we could optimize it more by finding the y_seg values automatically for you. 

piecewise_simulation = function(N = 1, log_f, y_seg, ...) {
  # assumes log_f is vectorized function and that f is log-concave
  # extra arguments to log_f can be given through ... 
  # also assumes that f is propto real target and f : I->(0,inf)
  # find log(f(x))' via numderiv:
  
  #need at least 2 given points, because we need two to make at least one z point. 
  stopifnot(length(y_seg) >= 2, all(diff(y_seg) > 0))
  
  y_seg = sort(y_seg)
  m = length(y_seg)
  a = vapply(y_seg,function(y) grad(log_f, x = y, ...),numeric(1)) #grad doesn't like vectorized stuff
  b = log_f(y_seg, ...) - a * y_seg
  z = diff(b)/(-diff(a))
  z = c(0, z, Inf)
  
  # technically this could be done vectorization mode, but cba.

  Q = numeric(m)
  for (i in 1:m) {
    Q[i] = exp(b[i]) * (exp(a[i] * z[i+1]) - exp(a[i] * z[i])) / a[i]
  }
  Q = cumsum(Q) #this has Q = [Q1, ..., d]
  
  u0 = as.numeric(Q[m]) * runif(N)
  u = runif(N)
  
  # the below is genuinely horrendously slow dogshit, should be vectorized. 
  I = matrix(NA, nrow=m, ncol = N)
  I[1, ] = (u0 < Q[1])
  if (m > 2) {
    for (i in 2:(m - 1)) {
      I[i, ] = (u0 >= Q[i - 1]) & (u0 < Q[i])
    }
  }
  I[m, ] = (u0 >= Q[m-1])
  
  x = numeric(N)
  accept = logical(N)
  
  x[I[1,]] = log(a[1] * exp(-b[1]) * u0[I[1,]] + exp(a[1] * z[1])) / a[1]
  accept[I[1,]] = log(u[I[1,]]) <= log_f(x[I[1,]], ...) - a[1] * x[I[1,]] - b[1]
  for (i in 2:length(y_seg)) {
    x[I[i,]] = log(a[i] * exp(-b[i]) * (u0[I[i,]] - Q[i-1]) + exp(a[i] * z[i])) / a[i]
    accept[I[i,]] = log(u[I[i,]]) <= log_f(x[I[i,]], ...) - a[i] * x[I[i,]] - b[i]
  }
  x[accept]
}







# Marie's attempt

# function from excercise, vectorized and takes x and z vectors from the data file poisson.csv
# this function could be a little faster using matrix rowsum stuff on the second term
# but that seems kinda unreadable ngl, so i like this. 
log_f = function(y, x_dat, z_dat) {
  stopifnot(all(is.integer(z_dat), na.rm = TRUE))
  stopifnot(all(x_dat > 0, na.rm = TRUE))
  stopifnot(all(y >= 0, na.rm = TRUE))
  # vectorized version so it's a little less annoying to deal with
  sxz = sum(x_dat*z_dat)
  tmp = function(y) {y*sxz - sum(exp(y*x_dat))}
  vals = vapply(y, tmp, numeric(1))
  vals
}

# analytical derivative of log_f
d_log_f <- function(y, x_dat, z_dat) {
  tmp = function(y) {sum(z_dat * x_dat) - sum(x_dat * exp(y * x_dat))}
  dvals = sapply(y, tmp)
  dvals
}


piecewise_simulation2 = function(N = 1, log_f, y_seg, ...) {
  # assumes log_f is vectorized function and that f is log-concave
  # extra arguments to log_f can be given through ... 
  # also assumes that f is propto real target and f : I->(0,inf)
  # find log(f(x))' via numderiv:
  
  #need at least 2 given points, because we need two to make at least one z point. 
  stopifnot(length(y_seg) >= 2, all(diff(y_seg) > 0)) 
  
  y_seg = sort(y_seg) # HERE WE DEFINE WHICH SEQUENCE OF Y TO INPUT
  m = length(y_seg) # NUMBER OF POINTS
  a = d_log_f(y_seq) # think faster to just use the defined function 
  b = log_f(y_seg, ...) - a * y_seg
  z = diff(b)/(-diff(a))
  z = c(0, z, Inf)
  
  Q = exp(b) * (exp(a * z[-1]) - exp(a * z[-length(z)])) / a
  Q = cumsum(Q)
  
  u0 = as.numeric(Q[m]) * runif(N)
  u = runif(N)
  
  I <- findInterval(u0, Q) + 1
  
  x = numeric(N)
  accept = logical(N)
  
  x = log(a[I] * exp(-b[I]) * (u0[I] - c(0, Q)[I]) + exp(a[I] * z[I])) / a[I]
  
  accept = log(u[I]) <= log_f(x[I])-a[I]*x[I]-b[I]
 
}

## Benchmarking

poisson <- read.csv("data/poisson.csv")
Z = poisson$z
X = poisson$x



