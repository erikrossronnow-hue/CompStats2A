# analytical derivative of log_f
d_log_f <- function(y, x_dat, z_dat) {
  tmp = function(y) {sum(z_dat * x_dat) - sum(x_dat * exp(y * x_dat))}
  dvals = sapply(y, tmp)
  dvals
}
