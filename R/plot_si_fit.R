#' Plot the epidemic curve and fitted serial interval distribution
#'
#' @param dat data frame containing a variable x with index case to case (ICC) intervals
#' @param mean mean of serial interval distribution
#' @param sd standard deviation of serial interval distribution
#' @param weights numeric vector of weights based on transmission route
#' @param dist string; assumed distribution of the serial interval; takes "normal" or "gamma";
#' defaults to "normal".
#' @param scaling_factor numeric; scales the density to better match the height of the histogram; defaults to 1.
#' @param x_scale numeric; scaling factor to control where the mean value is labelled along the x-axis. The mean value is positioned at `mean + (max(dat) * x_scale)`.
#' @return ggplot object
#' @export
#' @import ggplot2
#' @importFrom stats density
#'
plot_si_fit <- function(dat, mean, sd, weights, dist="normal", scaling_factor = 1,
                        x_scale = 0.04){

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
      theme_minimal() +
      # Add dashed vertical line at mean
      geom_vline(xintercept = mean, linetype = "dashed", color = "black") +
      # Annotate the mean value on the x-axis
      annotate("text", x = mean + (max(dat) * x_scale), y = min(density(dat)$y) - 0.001,
               label = paste0(round(mean, 1)), size = 3, color = "black", vjust = 1)
  }

  return(p)
}
