set.seed(1234)

# Parameters for simulation
N <- 500           # Number of observations
true_mean <- 15    # True mean serial interval (days)
true_sd <- 3       # True standard deviation (days)
route_weights <- c(0.2, 0.5, 0.2, 0.1)  # Weights for transmission routes

# Generate data for different transmission routes
CP <- fdrtool::rhalfnorm((route_weights[1]*N), theta=sqrt(pi/2)/(sqrt(2)*true_sd))  # Co-Primary
PS <- rnorm(route_weights[2]*N, mean=true_mean, sd=true_sd)                # Primary-Secondary
PT <- rnorm(route_weights[3]*N, mean=2*true_mean, sd=sqrt(2)*true_sd)      # Primary-Tertiary
PQ <- rnorm(route_weights[4]*N, mean=3*true_mean, sd=sqrt(3)*true_sd)      # Primary-Quaternary

# Combine and round to days
sim_icc_intervals <- round(c(CP, PS, PT, PQ))

