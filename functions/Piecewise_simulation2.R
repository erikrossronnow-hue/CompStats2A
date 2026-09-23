piecewise_simulation2 = function(N = 1, log_f, y_seq, ...) {
  # assumes log_f is vectorized function and that f is log-concave
  # extra arguments to log_f can be given through ... 
  # also assumes that f is propto real target and f : I->(0,inf)
  # find log(f(x))' via numderiv:
  
  #need at least 2 given points, because we need two to make at least one z point. 
  stopifnot(length(y_seq) >= 2, all(diff(y_seq) > 0)) 
  
  y_seq = sort(y_seq) # HERE WE DEFINE WHICH SEQUENCE OF Y TO INPUT
  m = length(y_seq) # NUMBER OF POINTS
  a = d_log_f(y_seq, X, Z) # think faster to just use the defined function 
  b = log_f(y_seq, X, Z) - a * y_seq
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
