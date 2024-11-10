## Functions to calculate exact Rt
## Code by Ed Baskerville
## June 13, 2020

#' Function to compute the probability, $p(u)$, that an individual is infectious u days after the moment of infection.
#'
#' @param sigma rate of leaving compartment E
#' @param gamma rate of leaving compartment I
#' @returns a function that computes the probability of infectiousnessu days after the moment of infection
#' @export
make_p_infectious <- function(sigma, gamma) {
  if(sigma == gamma) {
    rtn <- function(s) {
      gamma * s * exp(-gamma * s)
    }
  } else {
    rtn <- function(s) {
      sigma * (exp(-sigma * s) - exp(-gamma * s)) / (gamma - sigma)
    }
  }

  return(rtn)
}


#' Function to integrate SEIR
#'
#' @param beta_t Function that returns beta (t) for any time t.
#' @param sigma Rate of leaving latent class.
#' @param gamma Rate of leaving infected class.
#' @param N Total population size.
#' @param T_final Final time step.
#' @param E0 Initial number of individuals in exposed state.
#' @param I0 Initial number of individuals in infectious state.
#' @param dt Size of time step over which to integrate SEIR; defaults to `dt = 0.01`.
#' @returns Data frame with number of individuals in each state for each time point.
#' @export
#' @importFrom deSolve ode
integrate_seir <- function(beta_t,
                           sigma,
                           gamma,
                           N,
                           T_final,
                           E0,
                           I0,
                           dt = 0.01) {
  S0 <- N - E0 - I0
  ddt_ode <- function(t, y, params) {
    S <- y[1]
    E <- y[2]
    I <- y[3]

    dSE <- beta_t(t) * S * I / N
    dEI <- sigma * E
    dIR <- gamma * I

    list(c(S = -dSE, E = dSE - dEI, I = dEI - dIR))
  }

  run_ode <- as.data.frame(ode(
    c(
      S = S0, E = E0, I = I0
    ), seq(0, T_final, dt), ddt_ode, list(), method = 'rk4'
  ))
  return(run_ode)
}


# # Example
# arnaught <- 2
# t_E <- 2
# t_I <- 4
# beta_t <- function(t) { arnaught / t_I } # Can replace this with a beta(t) with changes
# sigma <- 1 / t_E
# gamma <- 1 / t_I
# N <- 1
# T <- 80
# E0 <- 0.001
# I0 <- 0
# ode_output <- as.data.frame(integrate_seir(beta_t, sigma, gamma, N, T, E0, 0))
# ggplot(ode_output, aes(x = time, y = S / N)) + geom_line()



#' Calculate the instantaneous case reproductive number using the deterministic SEIR output
#'
#' @param beta_t Function that returns beta (t) for any time t.
#' @param sigma Rate of leaving latent class.
#' @param gamma Rate of leaving infected class.
#' @param N Total population size.
#' @param T_final Final time.
#' @param E0 Initial number of individuals in exposed state.
#' @param I0 Initial number of individuals in infectious state.
#' @param dt_int Integration time step; defaults to `dt_int = 0.1`.
#' @param dt_Rt Case reproduction number output time step; defaults to `dt_Rt = 1`.
#' @param T_inf Final time in SEIR integration. Should be at least one generation interval higher than T_final.
#' @returns Named list with two data frames:
#'  1. `case_rt_df` with estimated case reproduction number for every time point, and
#'  2. `ode_output` with the output from the compartmental model used to estimate Rt
#' @export
#' @importFrom stats approxfun
integrate_Rt <- function(beta_t,
                         sigma,
                         gamma,
                         N,
                         T_final,
                         E0,
                         I0,
                         dt_int = 0.1,
                         dt_Rt = 1,
                         T_inf = 1000) {

  S0 <- (N - E0 - I0)
  T_inf = T_final + 100
  if(any(is.na(beta_t(0:T_inf)))){
    stop('beta_t must return values from 0:(T+100). Either decrease T or increase the range of beta_t.')
    }

  ode_output <- as.data.frame(
    integrate_seir(beta_t = beta_t,
                   sigma = sigma,
                   gamma = gamma,
                   N = N,
                   T_final = T_inf,
                   E0 = E0,
                   I0 = I0)
  )

  SoverN <- approxfun(ode_output$time, ode_output$S / N)

  p_infectious <- make_p_infectious(sigma, gamma)

  times <- seq(0, T_final, dt_Rt)
  # Rt <- sapply(times, function(t) {
  #   integrate(function(u) {
  #     beta_t(u) * p_infectious(u - t) * SoverN(u)
  #   }, t, T_inf)$value
  # })

  Rt <- sapply(times, function(time) {
    int_func <- function(u) {
      beta_val <- beta_t(u)
      p_inf_val <- p_infectious(u - time)
      s_n_val <- SoverN(u)
      # cat("u:", u, "\n",
      #     "beta_t:", beta_val, "\n",
      #     "p_infectious:", p_inf_val, "\n",
      #     "SoverN:", s_n_val, "\n")
      beta_val * p_inf_val * s_n_val
    }

    integrate(int_func, time, T_inf)$value
  })

  rtn <- list(
    case_rt_df = data.frame(time = times, R_case = Rt),
    ode_output = ode_output
  )

  return(rtn)
}
