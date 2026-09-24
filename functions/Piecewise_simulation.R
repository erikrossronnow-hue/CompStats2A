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
  
  # Find a, b and boundary points z 
  a = vapply(y_seg,function(y) grad(log_f, x = y, ...),numeric(1)) 
  b = log_f(y_seg, ...) - a * y_seg
  z = diff(b)/(-diff(a))
  z = c(0, z, Inf)
  
  Q = numeric(m)
  for (i in 1:m) {
    Q[i] = exp(b[i]) * (exp(a[i] * z[i+1]) - exp(a[i] * z[i])) / a[i]
  }
  Q = cumsum(Q)
  
  u0 = as.numeric(Q[m]) * runif(N)
  u = runif(N)
  
  # find interval
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
  
  # Sample x and Accept/Reject
  x[I[1,]] = log(a[1] * exp(-b[1]) * u0[I[1,]] + exp(a[1] * z[1])) / a[1]
  accept[I[1,]] = log(u[I[1,]]) <= log_f(x[I[1,]], ...) - a[1] * x[I[1,]] - b[1]
  
  for (i in 2:length(y_seg)) {
    x[I[i,]] = log(a[i] * exp(-b[i]) * (u0[I[i,]] - Q[i-1]) + exp(a[i] * z[i])) / a[i]
    accept[I[i,]] = log(u[I[i,]]) <= log_f(x[I[i,]], ...) - a[i] * x[I[i,]] - b[i]
  }
  
  x[accept]
}
