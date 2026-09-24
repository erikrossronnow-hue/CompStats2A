piecewise_simulation2 = function(N = 1, log_f, y_seg, ...) {
  # assumes log_f is vectorized function and that f is log-concave
  # extra arguments to log_f can be given through ... 
  # also assumes that f is propto real target and f : I->(0,inf)
  
  #need at least 2 given points, because we need two to make at least one z point. 
  stopifnot(length(y_seg) >= 2, all(diff(y_seg) > 0)) 
  
  y_seg = sort(y_seg) 
  m = length(y_seg)
  
  # Find a, b and boundary points z 
  a = d_log_f(y_seg, ...) 
  b = log_f(y_seg, ...) - a * y_seg
  z = diff(b)/(-diff(a))
  z = c(0, z, Inf)
  
  Qi = exp(b) * (exp(a * z[-1]) - exp(a * z[-length(z)])) / a
  Q = cumsum(Qi)
  Q0 = c(0, Q[-m]) 
  
  u0 = Q[m] * runif(N)
  u = runif(N)
  
  # Intervals
  I = findInterval(u0, Q) + 1
  I[I > m] = m
  
  # Sample x
  x = log(a[I] * exp(-b[I]) * (u0 - Q0[I]) + exp(a[I] * z[I])) / a[I]
  
  # Accept/Reject
  accept = log(u) <= log_f(x, ...) - a[I] * x - b[I]
  
  x[accept]
  
}
