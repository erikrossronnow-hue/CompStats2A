piecewise_simulation2 = function(N = 1, log_f, y_seg, ...) {
  # assumes log_f is vectorized function and that f is log-concave
  # extra arguments to log_f can be given through ... 
  # also assumes that f is propto real target and f : I->(0,inf)
  # find log(f(x))' via numderiv:
  
  #need at least 2 given points, because we need two to make at least one z point. 
  stopifnot(length(y_seg) >= 2, all(diff(y_seg) > 0)) 
  
  y_seg = sort(y_seg) # HERE WE DEFINE WHICH SEQUENCE OF Y TO INPUT
  m = length(y_seg) # NUMBER OF POINTS
  a = d_log_f(y_seg, ...) # think faster to just use the defined function 
  b = log_f(y_seg, ...) - a * y_seg
  z = diff(b)/(-diff(a))
  z = c(0, z, Inf)
  
  Qi = exp(b) * (exp(a * z[-1]) - exp(a * z[-length(z)])) / a
  Q = cumsum(Qi)
  Q0 = c(0, Q[-m]) 
  
  u0 = Q[m] * runif(N)
  u = runif(N)
  
  I = findInterval(u0, Q) + 1
  I[I > m] = m
  
  x = log(a[I] * exp(-b[I]) * (u0 - Q0[I]) + exp(a[I] * z[I])) / a[I]
  
  accept = log(u) <= log_f(x, ...) - a[I] * x - b[I]
  
  x[accept]
  
}
