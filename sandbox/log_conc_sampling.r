# comp stat rejection sampling sandbox script
rm(list = ls())
setwd("~/Desktop/CompStat2A")
getwd()
library(data.table)
library(numDeriv)
library(profvis)
library(bench)
# if u have downloaded this script u need to grab the poisson.csv file urself!
#poisson = read.csv("poisson.csv")
poisson = read.csv("data/poisson.csv")
DT = as.data.table(poisson)
DT[, z := as.integer(z)]

#load the necessary functions, not many of them so no fancy setup. 

#source("functions\log_f.r")
#source("functions\Piecewise_simulation.R")

source("log_f.r")
source("d_log_f.r")
source("Piecewise_simulation.R")
source("Piecewise_simulation2.R")

# want to plot the function so we can see how many breakpoints we need approx:
# chatmaxxing btw to make it look niceeeee :))
yy = seq(0.01, 0.5, by = 0.001)
log_vals <- log_f(yy, x_dat = DT$x, z_dat = DT$z)
plot(
  yy,
  exp(log_vals - max(log_vals)),
  type = "l",
  main = "Rescaled target density",
  xlab = expression(y),
  ylab = expression(tilde(f)(y) / max(tilde(f)))
)

yy <- seq(0.01, 0.5, by = 0.001)

log_vals <- log_f(yy, x_dat = DT$x, z_dat = DT$z)

plot(
  yy,
  log_vals,
  type = "l",
  main = "Log target density",
  xlab = expression(y),
  ylab = expression(log(f(y)))
)
# unironically just a bell curve thingie, completely fine to just have 
# m=2 here, one segments from y in (0,0.25) and another from(0.25, inf)
# so chose y_1 = 0.18 and y_2 = 0.3 for each side of the curve.

# trial 1
set.seed(1)
y_seg = c(0.18, 0.3) # note that these y's are on each side of the bell curve 
sim = piecewise_simulation(N=1000, log_f, y_seg=y_seg, x_dat = DT$x, z_dat = DT$z)
hist(sim, breaks=50)


# trial 2 
sim2 = piecewise_simulation2(N=1000, log_f, y_seg=y_seg, x_dat = DT$x, z_dat = DT$z)
hist(sim2, breaks = 50)

# benchmarking: 
N_bench = 10000
y_seg = c(0.18, 0.30)

bm_piecewise = bench::mark(
  
  piecewise = {
    set.seed(420)
    
    piecewise_simulation(
      N = N_bench,
      log_f = log_f,
      y_seg = y_seg,
      x_dat = DT$x,
      z_dat = DT$z
    )
  },
  
  piecewise2 = {
    piecewise_simulation2(
      N=N_bench, 
      log_f, 
      y_seg=y_seg, 
      x_dat = DT$x, 
      z_dat = DT$z)
  },
  
  iterations = 100,
  check = FALSE,
  memory = TRUE
)

bm_piecewise
plot(bm_piecewise)


set.seed(123) 
#acceptance rate stuff
y = piecewise_simulation(
  N = N_bench,
  log_f = log_f,
  y_seg = c(0.18, 0.30),
  x_dat = DT$x,
  z_dat = DT$z
)

y2 = piecewise_simulation2(
  N = N_bench,
  log_f = log_f,
  y_seg = c(0.18, 0.30),
  x_dat = DT$x,
  z_dat = DT$z
)

length(y) / N_bench # around 0.7 
length(y2) / N_bench


# profiling

p1 = profvis::profvis({
  
  set.seed(123)
  
  y <- piecewise_simulation(
    N = 5e6,
    log_f = log_f,
    y_seg = c(0.18, 0.30),
    x_dat = DT$x,
    z_dat = DT$z
  )
  
}, interval = 0.001)


# more or less confirms that the stuff we knew was insane dogshit, was in fact
# insane dogshit. 
p1


p2 = profvis::profvis({
  
  set.seed(123)
  
  y <- piecewise_simulation(
    N = 5e6,
    log_f = log_f,
    y_seg = c(0.1, 0.18, 0.30, 0.4), #finer grid, does it dominate profile?
    x_dat = DT$x,
    z_dat = DT$z
  )
  
}, interval = 0.001)
p2

p3 = profvis::profvis({
  
  set.seed(123)

  y2 = piecewise_simulation2(
    N = 5e6,
    log_f = log_f,
    y_seg = c(0.18, 0.30),
    x_dat = DT$x,
    z_dat = DT$z
  )
  
  
}, interval = 0.001)



# mega chat stuff in terms of looking at acceptance rates versus m:
library(bench)
library(data.table)
library(ggplot2)

N_bench <- 100000
m_values <- c(2, 3, 4, 5, 6, 8, 10)

# These are only convenient initial segment placements.
# Replace them with your own optimal/adaptive y_seg values if desired.
y_seg_list <- lapply(m_values, function(m) {
  seq(0.10, 0.40, length.out = m)
})

benchmark_results <- rbindlist(lapply(seq_along(m_values), function(j) {
  
  m <- m_values[j]
  y_seg <- y_seg_list[[j]]
  
  bm <- bench::mark(
    piecewise = {
      set.seed(420)
      
      piecewise_simulation(
        N = N_bench,
        log_f = log_f,
        y_seg = y_seg,
        x_dat = DT$x,
        z_dat = DT$z
      )
    },
    iterations = 30,
    check = FALSE,
    memory = TRUE
  )
  
  # Separate deterministic run to obtain the acceptance rate
  set.seed(420)
  
  y_accepted <- piecewise_simulation(
    N = N_bench,
    log_f = log_f,
    y_seg = y_seg,
    x_dat = DT$x,
    z_dat = DT$z
  )
  
  data.table(
    m = m,
    median_ms = as.numeric(bm$median) * 1000,
    min_ms = as.numeric(bm$min) * 1000,
    max_ms = as.numeric(max(bm$time[[1]])) * 1000,
    memory_mb = as.numeric(bm$mem_alloc) / 1024^2,
    acceptance_rate = length(y_accepted) / N_bench
  )
}))

benchmark_results


ggplot(benchmark_results, aes(x = m, y = median_ms)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = min_ms, ymax = max_ms),
    width = 0.15,
    alpha = 0.6
  ) +
  scale_x_continuous(breaks = m_values) +
  labs(
    title = "Piecewise rejection sampler runtime by number of segments",
    x = "Number of envelope segments, m",
    y = "Runtime per 100,000 proposals (ms)"
  ) +
  theme_minimal(base_size = 13)

ggplot(benchmark_results, aes(x = m, y = acceptance_rate)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = m_values) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  labs(
    title = "Acceptance rate by number of envelope segments",
    x = "Number of envelope segments, m",
    y = "Acceptance rate"
  ) +
  theme_minimal(base_size = 13)


plot_data <- melt(
  benchmark_results,
  id.vars = "m",
  measure.vars = c("median_ms", "acceptance_rate"),
  variable.name = "measure",
  value.name = "value"
)

ggplot(plot_data, aes(x = m, y = value)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ measure, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = m_values) +
  labs(
    title = "Effect of envelope complexity on sampler performance",
    x = "Number of envelope segments, m",
    y = NULL
  ) +
  theme_minimal(base_size = 13)

benchmark_results[, ms_per_accepted_draw :=
                    median_ms / (N_bench * acceptance_rate)
]



# Grid for plotting
yy <- seq(0.001, 0.5, by = 0.001)

# True log-density
log_vals <- log_f(
  yy,
  x_dat = DT$x,
  z_dat = DT$z
)

# Rescale true density
f_true <- exp(log_vals - max(log_vals))


yy <- seq(0.01, 0.5, by = 0.001)

log_vals <- log_f(yy, x_dat = DT$x, z_dat = DT$z)

y_seg <- c(0.18, 0.30)

# Tangent slopes
a <- sapply(
  y_seg,
  d_log_f,
  x_dat = DT$x,
  z_dat = DT$z
)

# Tangent intercepts
b <- log_f(
  y_seg,
  x_dat = DT$x,
  z_dat = DT$z
) - a * y_seg

# Intersection of the two tangent lines
z <- (b[2] - b[1]) / (a[1] - a[2])

z

left <- yy <= z
right <- yy >= z

tangent_left <- a[1] * yy + b[1]
tangent_right <- a[2] * yy + b[2]

plot(
  yy,
  log_vals,
  type = "l",
  lwd = 2,
  xlab = expression(y),
  ylab = expression(log(f(y))),
  main = "Log-density with tangent envelope"
)

lines(
  yy[left],
  tangent_left[left],
  col = "red",
  lwd = 2
)

lines(
  yy[right],
  tangent_right[right],
  col = "red",
  lwd = 2
)

abline(v = z, lty = 2)

