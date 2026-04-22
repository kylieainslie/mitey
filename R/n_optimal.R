n_optimal1 <- function(dat) {
  as.integer(quantile(dat, 0.95) / mean(dat))
}


n_optimal2 <- function(dat) {
  as.integer(sd(dat) / mean(dat)) + 2L  # +2 car minimum 2 routes
}

n_optimal3 <- function(dat) {
  as.integer(max(dat) / sd(dat))
}

n_optimal1(icc_Influenza_Canada)
n_optimal2(icc_Influenza_Canada)
n_optimal3(icc_Influenza_Canada)
