setwd(r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\mitey\mitey\R)")

source("weighted_var.R")
source("wt_loglik.R")

source("flower.R")
source("fupper.R")
source("f0.R")
source("integrate_component.R")

source("f_norm.R")
source("f_gam.R")
source("f_lnorm.R")     
source("si_estim.R")
source("plot_si_fit.R")

# smoke test
set.seed(123)

dat_test <- c(
  0, 1, 1, 2, 2, 3, 4, 5, 6, 8,
  10, 12, 14, 16, 18
)

fit_lnorm <- si_estim(
  dat = dat_test,
  n = 3,
  dist = "lognormal"
)

fit_lnorm

# test2

set.seed(2026)

mean_si <- 6
sd_si <- 3

sdlog <- sqrt(log(1 + sd_si^2 / mean_si^2))
meanlog <- log(mean_si) - 0.5 * sdlog^2

n_cp <- 10
n_ps <- 25
n_pt <- 10
n_pq <- 5

cp <- abs(
  rlnorm(n_cp, meanlog, sdlog) -
    rlnorm(n_cp, meanlog, sdlog)
)

ps <- rlnorm(n_ps, meanlog, sdlog)

pt <- rlnorm(n_pt, meanlog, sdlog) +
  rlnorm(n_pt, meanlog, sdlog)

pq <- rlnorm(n_pq, meanlog, sdlog) +
  rlnorm(n_pq, meanlog, sdlog) +
  rlnorm(n_pq, meanlog, sdlog)

dat_lnorm <- round(c(cp, ps, pt, pq))
dat_lnorm <- pmax(dat_lnorm, 0)

fit_lnorm <- si_estim(
  dat = dat_lnorm,
  n = 5,
  dist = "lognormal"
)
fit_lnorm

library(ggplot2)
p_lnorm <- plot_si_fit_result(
  si_result = fit_lnorm,
  dat = dat_lnorm,
  dist = "lognormal",
  scaling_factor = 1
)

print(p_lnorm)
ggplot2::ggsave(
  filename = "test_lognormal_fit.png",
  plot = p_lnorm,
  width = 7,
  height = 5,
  dpi = 300
)

# test original functions
library(fdrtool)
fit_norm <- si_estim(
  dat = dat_test,
  n = 3,
  dist = "normal"
)

fit_gam <- si_estim(
  dat = dat_test,
  n = 3,
  dist = "gamma"
)

stopifnot(is.finite(fit_norm$mean))
stopifnot(is.finite(fit_norm$sd))
stopifnot(length(fit_norm$wts) == 7)

stopifnot(is.finite(fit_gam$mean))
stopifnot(is.finite(fit_gam$sd))
stopifnot(length(fit_gam$wts) == 4)
