# Rt vailidation with EpiLPS

# Load packages
library(devtools)
load_all() # load mitey
#devtools::install_github("oswaldogressani/EpiLPS")
library(EpiLPS) # load Epi

# set seed
set.seed(1234)

# specify serial interval distribution
mean_si <- 8
sd_si <- 5.66

si_spec <- Idist(mean = mean_si, sd = sd_si, dist = "gamma")
plot(si_spec)

# generate data
t_end <- 40
datasim <- episim(si = si_spec$pvec, Rpattern = 5, endepi = t_end,
                  dist = "negbin", overdisp = 15, plotsim = TRUE)

# Estimate Rt w EpiLPS
# fit <- estimR(incidence = datasim$y, si = si_spec$pvec)
# fit

fitmcmc <- estimRmcmc(incidence = datasim$y, si = si_spec$pvec,
                      CoriR = TRUE, WTR = TRUE,
                      niter = 5000, burnin = 2000)
summary(fitmcmc)

# Estimate Rt with mitey::rt_estim()

t_vec <- seq(1, length(datasim$y), by = 1)
inc_dat <- data.frame(onset_date = t_vec, inc = datasim$y)

rt_estimated <- rt_estim_w_boot(
  inc_dat = inc_dat,
  mean_si = mean_si,
  sd_si = sd_si,
  dist_si = "gamma",
  n_bootstrap = 100
)

rt_WL <- ifelse(is.nan(rt_estimated$rt), NA, rt_estimated$rt)
rt_WL <- ifelse(is.infinite(rt_estimated$rt) | rt_estimated$rt > 15, NA, rt_estimated$rt)
rt_WL_smooth <- rollmean(rt_WL, k = 7, fill = NA, na.rm = TRUE)
rt_WL_shift <- c(rep(NA, mean_si), rt_WL_smooth)

# Plot
tt <- seq(mean_si + 1, t_end, by = 1)
Rtrue <- sapply(tt, datasim$Rtrue)
plot(tt, Rtrue, type = "l", xlab = "Time", ylab = "R", ylim = c(0,4),
     lwd = 2)
lines(tt, fitmcmc$RLPS$R[-(1:mean_si)], col = "red", lwd = 2)
lines(tt, fitmcmc$RCori$`Mean(R)`, col = "blue", lwd = 2)
lines(tt, fitmcmc$RWT$`Mean(R)`[-1], col = "purple", lwd = 2)
lines(tt, rt_WL_smooth[tt], col = "green", lwd = 2)
lines(tt, rt_WL_shift[tt], col = "green", lty = 2, lwd = 2)
legend("topright", col = c("black","red","blue", "purple", "green", "green"),
       c("True R","EpiLPS","EpiEstim", "WT", "WL", "WL (shifted)"),
       bty = "n", lty = c(1,1,1,1,1,2))


