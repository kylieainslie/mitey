//
// This Stan program defines a model to correct for right truncation when using the Wallinga and Teunis method to estimate case reproduction number (Rt). We assume that the generation time is unknown, so use the serial interval to approximate the generation time. We assume a normally distributed serial interval. We aim to estimate the expected Rt when correcting for secondary cases that have not yet been observed.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  real right_truncation_correction(real t, real T, real mean_serial_interval, real sd_serial_interval) {
    // Probability that the serial interval is less than the time lag (T - t)
    real prob = Phi((T - t - mean_serial_interval) / sd_serial_interval);
    return 1.0 / prob; // Correction factor as per Cauchemez et al.
  }
}

// The input data.
data {
  int<lower=0> N;                    // Number of observations
  vector[N] time_since_infection;    // Time since each primary case was observed
  vector[N] Rt_obs;                  // Observed Rt values (can be derived from Wallinga-Teunis method)
  real<lower=0> T;                   // Current time point (most recent observation)
  real<lower=0> mean_serial_interval;// Mean serial interval (used as a proxy for generation interval)
  real<lower=0> sd_serial_interval;  // Standard deviation of the serial interval
}

// The parameters accepted by the model. Our model
// accepts one parameter 'Rt'.
parameters {
  real<lower=0> Rt[N];   // True Rt values
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  Rt ~ normal(Rt_obs, 0.1); // Prior for Rt

  // Likelihood function incorporating right truncation as per Cauchemez et al.
  for (i in 1:N) {
    // Correction factor using the serial interval as a proxy for the generation interval
    real truncation_correction = right_truncation_correction(time_since_infection[i], T, mean_serial_interval, sd_serial_interval);

    // Adjusted likelihood incorporating right truncation
    target += normal_lpdf(Rt[i] | Rt_obs[i] * truncation_correction, 0.1);
  }
}
