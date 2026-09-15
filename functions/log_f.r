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
