#' Simulate a deterministic ODE approximation of a continuous-time, discrete-state stochastic S(E)IR model.
#'
#' Original code created by Ed Baskerville (15 April 2020)
#'
#' This model has no age structure.
#'
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0. Average number of new infections produced by an infection in a susceptible population. A scalar or a vector of length `n_t + 1`, which specifies R0 at the start of each timestep. R0 is linearly interpolated between timesteps.
#' @param t_E Mean latent period. If set to 0, the model reduces to an SIR.
#' @param t_I Mean duration of infectiousness.
#' @param N integer; Total population size.
#' @param S_init integer; number of individuals starting in the susceptible state.
#' @param E_init integer; number of individuals starting in the exposed state.
#' @param I_init interger; number of individuals starting in the infectious state.
#' @param n_t Number of units of time to simulate.
#' @param n_steps_per_t integer; number of time steps, defaults to 1.
#' @returns data frame with number of individuals transitioning between each state per time step
#' @export
#' @import deSolve
#' @import tidyr
#'
simulate_seir_ode <- function(
    arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t,
    n_steps_per_t = 1 # Ignored; included so the function signature matches stochastic version
) {

  check_args(
    arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t, n_steps_per_t
  )

  beta <- construct_beta(arnaught, t_I, n_t)
  d_dt <- function(t, y, params) {
    dS <- y['S'] * beta(t) * y['I'] / N
    dIR <- y['I'] / t_I

    if(t_E > 0) {
      # SEIR model
      dEI <- y['E'] / t_E
      list(c(
        S = -dS,
        E = dS - dEI,
        I = dEI - dIR,
        R = dIR,
        cum_dS = dS,
        cum_dEI = dEI
      ), NULL)
    }
    else {
      # SIR model
      list(c(
        S = -dS,
        E = 0,
        I = dS - dIR,
        R = dIR,
        cum_dS = dS,
        cum_dEI = dS
      ), NULL)
    }
  }

  y_init <- c(
    S = S_init,
    E = if(t_E > 0) E_init else 0,
    I = if(t_E > 0) I_init else E_init + I_init,
    R = 0,
    cum_dS = 0,
    cum_dEI = 0
  )

  # automatic ode solver is lsoda, an "automatic stiff/non-stiff solver"
  results_ode <- ode(y = y_init, times = 0:n_t, func = d_dt, parms = NULL)

  rtn <- as.data.frame(results_ode) %>%
    mutate(dS = cum_dS - lag(cum_dS, 1)) %>%
    mutate(dEI = cum_dEI - lag(cum_dEI, 1)) %>%
    mutate(dIR = R - lag(R, 1))

  return(rtn)
}
