# Helper functions for generating simulated data from a S(E)IR model

#' Function to check if input number is an integer
#' @param x numeric; value to check whether it is a whole number
#' @returns logical expression; returns TRUE if x is a whole number
#' @export
is.wholenumber <- function(x) {
  if(is.character(x) || is.logical(x)){
    stop("x must be numerical")
  }
  rtn <- x == round(x)
  return(rtn)
}


#' Function to check the arguments of `simulate_seir_ode()`.
#'
#' The function checks that the initial values for the different states are whole numbers and that $R_0$, $t_E$ and $t_I$ are numeric.
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0. Average number of new infections produced by an infection in a susceptible population. A scalar or a vector of length `n_t + 1`, which specifies R0 at the start of each timestep. R0 is linearly interpolated between timesteps.
#' @param t_E numeric; Mean latent period. If set to 0, the model reduces to an SIR.
#' @param t_I numeric; Mean duration of infectiousness.
#' @param N integer; Total population size.
#' @param S_init integer; number of individuals starting in the susceptible state.
#' @param E_init integer; number of individuals starting in the exposed state.
#' @param I_init interger; number of individuals starting in the infectious state.
#' @param n_t integer; Number of units of time to simulate.
#' @param n_steps_per_t integer; number of time steps, defaults to 1.
#' @returns error message
#' @export
check_args <- function(
    arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t, n_steps_per_t
) {
  # Check all input parameters
  stopifnot(is.wholenumber(N) && length(N) == 1 && N >= 1)
  stopifnot(is.wholenumber(n_t) && length(n_t) == 1 && n_t >= 1)
  stopifnot(
    is.wholenumber(n_steps_per_t) && length(n_steps_per_t) == 1 &&
      n_steps_per_t >= 1
  )
  stopifnot(
    is.numeric(arnaught) && arnaught > 0 &&
      (length(arnaught) == 1 || length(arnaught) == n_t + 1)
  )
  stopifnot(
    is.numeric(t_E) && length(t_E) == 1 && t_E >= 0
  )
  stopifnot(
    is.numeric(t_I) && length(t_I) == 1 && t_I > 0
  )
}

#' Function to construct the transmission rate from $R_0$, $t_I$ (mean infectious period), and $n_t$.
#'
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0. Average number of new infections produced by an infection in a susceptible population. A scalar or a vector of length `n_t + 1`, which specifies R0 at the start of each timestep. R0 is linearly interpolated between timesteps.
#' @param t_I numeric; Mean duration of infectiousness.
#' @param n_t integer; Number of units of time to simulate.
#' @returns a function that determines the transmission rate at each time point t
#' @export
#' @importFrom stats approxfun
construct_beta <- function(arnaught, t_I, n_t) {
  beta_t_all <- arnaught / t_I
  if(length(arnaught) == 1) {
    function(t) beta_t_all
  } else {
    approxfun(0:n_t, beta_t_all)
  }
}
