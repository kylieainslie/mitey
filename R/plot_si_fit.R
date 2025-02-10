#' Plot the epidemic curve and fitted serial interval distribution
#'
#' @param dat data frame containing a variable x with index case to case (ICC) intervals
#' @param mean mean of serial interval distribution
#' @param sd standard deviation of serial interval distribution
#' @param weights numeric vector of weights based on transmission route
#' @param dist string; assumed distribution of the serial interval; takes "normal" or "gamma";
#' defaults to "normal".
#' @param scaling_factor numeric; scales the density to better match the height of the histogram; defaults to 1.
#' @return ggplot object
#' @export
#' @import ggplot2
#'
plot_si_fit <- function(dat, mean, sd, weights, dist="normal", scaling_factor = 1){

  if (dist == "gamma"){

  # Define breaks
  breaks <- seq(min(dat) - 0.51, max(dat) + 0.51, by = 1)

  # Create the plot
  p <- ggplot(data = data.frame(x = dat), aes(x = .data$x)) +
    geom_histogram(aes(y = after_stat(.data$density) * scaling_factor), breaks = breaks, fill = "lightblue", color = "black") +
    stat_function(fun = f_gam,
                  args = list(w1 = weights[1], w2 = weights[2], w3 = weights[3],
                              mu = mean, sigma = sd),
                  color = "red", linetype = "dashed", size = 1) +
    labs(x = "Index-case to case interval (days)", y = "Probability") +
    theme_minimal()
  }

  if (dist == "normal"){
    # Define breaks for histogram
    breaks <- seq(min(dat) - 0.51, max(dat) + 0.51, by = 1)

    # Create the plot
    p <- ggplot(data = data.frame(x = dat), aes(x = .data$x)) +
      geom_histogram(aes(y = after_stat(.data$density) * scaling_factor), breaks = breaks, fill = "lightblue", color = "black") +
      stat_function(fun = f_norm, args = list(w1 = weights[1], w2 = weights[2], w3 = weights[3],
                                          mu = mean, sigma = sd),
                    color = "red", linetype = "solid", linewidth = 1) +
      labs(x = "Index-case to case interval (days)", y = "Density") +
      theme_minimal()
  }

  return(p)
}
