# Using a graph from Mellanby 1944 (Fig. 5, page 204) we will reconstruct the
# distribution of the number of mites during a scabies infection.
# We will then assess whether the right tail is better described by exponential
# decline or power law decline using the poweRlaw package.

# citation for poweRlaw package
# Colin S. Gillespie (2015). Fitting Heavy Tailed Distributions: The poweRlaw 
# Package. Journal of Statistical Software, 64(2), 1-16. 
# http://www.jstatsoft.org/v64/i02/.

# load required packages
library(poweRlaw)
library(tidyverse)

# create data set using Fig. 5 from Mellnaby 1944
mite_dist_df <- data.frame(
  days = c(1, 36, 50, 90, 100, 115, 125, 150, 200, 220, 240, 250, 260, 275, 290, 300),
  parasite_rate = c(1, 5, 50, 125, 150, 125, 100, 50, 10, 7, 10, 8, 7, 10, 8, 9)
)

# plot to check if it matches the figure
plot(parasite_rate ~ days, data = mite_dist_df)

# fit a power law distribution to the data
m1 = displ$new(mite_dist_df$parasite_rate)
m1$setPars(estimate_pars(m1))

# now fit exponential distribution
m2 = disexp$new(mite_dist_df$parasite_rate)
m2$setPars(estimate_pars(m2))

# plot both distributions to visually compare
plot(m1, ylab = "CDF")
lines(m1)
lines(m2, col = 2, lty = 2)

# it looks like exponential is a better fit, but let's compare using Vuong’s
# test statistic. We're performing a one-sided test to determine if the exponential
# distribution (m2) is a better fit than the power law (m1). NOTE: order matters
# in compare_distribution() for a one-sided test
compare_distributions(m2, m1)$p_one_sided
# [1] 0.05417215
# according to these results, an exponential dist is better a better fit to the 
# data.

# plot fitted distribution against data
plot(parasite_rate ~ days, data = mite_dist_df)

# fit data using loess function
rng <- range(mite_dist_df$days)
xx <- seq(rng[1], rng[2], length = 100)
lo <- loess(parasite_rate ~ days, mite_dist_df, span = 0.35)
lines(predict(lo, data.frame(days = xx)) ~ xx)

# fit data with gam
library(mgcv) # mgcv comes with R.  No need to install. Just load.
ga <- gam(parasite_rate ~ s(days, bs = "cs"), data = mite_dist_df)
lines(predict(ga, data.frame(days = xx)) ~ xx, col = 2, lty = 2)

legend("topleft", legend = c("loess", "gam"), lty = 1:2, col = 1:2)
