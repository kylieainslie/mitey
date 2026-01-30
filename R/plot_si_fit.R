#' Visualize Serial Interval Distribution Fit to Outbreak Data
#'
#' Creates a diagnostic plot showing the fitted serial interval mixture distribution
#' overlaid on a histogram of observed index case-to-case (ICC) intervals from outbreak
#' data.
#'
#' The function displays:
#' \itemize{
#'   \item \strong{Histogram}: Observed ICC intervals binned by day, representing the
#'         empirical distribution of time differences between symptom onset in the
#'         index case and all other cases in the outbreak
#'   \item \strong{Fitted curve}: The estimated mixture distribution combining different
#'         transmission routes (co-primary, primary-secondary, primary-tertiary, and
#'         primary-quaternary), weighted according to their estimated probabilities
#'   \item \strong{Reference line}: For normal distributions, a dashed vertical line
#'         indicates the estimated mean serial interval
#' }
#'
#'
#' @param dat numeric vector; the index case-to-case (ICC) intervals in days.
#'            These represent the time differences between symptom onset in the
#'            index case (case with earliest symptom onset) and each other case
#'            in the outbreak
#' @param mean numeric; the estimated mean of the serial interval distribution in days,
#'             typically obtained from \code{si_estim()}
#' @param sd numeric; the estimated standard deviation of the serial interval
#'           distribution in days, typically obtained from \code{si_estim()}
#' @param weights numeric vector; the estimated weights for different transmission
#'                route components. Length and interpretation depends on distribution:
#' \itemize{
#'   \item \strong{Normal distribution}: 4 weights corresponding to aggregated
#'         transmission routes (co-primary, primary-secondary, primary-tertiary,
#'         primary-quaternary)
#'   \item \strong{Gamma distribution}: 3 weights for the reduced component set
#' }
#' @param dist character; the distribution family used for serial interval estimation.
#'             Must be either "normal" (default) or "gamma". Should match the
#'             distribution used in the original \code{si_estim()} call
#' @param scaling_factor numeric; multiplicative factor to adjust the height of the
#'                       fitted density curve relative to the histogram. Values > 1
#'                       make the curve higher, values < 1 make it lower. Defaults to 1.
#'                       Useful when histogram and density have different scales.
#'
#' @return A \code{ggplot2} object that can be further customized or displayed.
#'         The plot includes appropriate axis labels, legend, and styling for
#'         publication-quality figures
#'
#' @seealso \code{\link{si_estim}} for serial interval estimation,
#'          \code{\link{f_norm}} and \code{\link{f_gam}} for the underlying
#'          mixture distribution functions
#'
#' @references
#' Vink MA, Bootsma MCJ, Wallinga J (2014). Serial intervals of respiratory infectious
#' diseases: A systematic review and analysis. American Journal of Epidemiology,
#' 180(9), 865-875.
#'
#' @export
#' @import ggplot2
#' @importFrom stats density
#'
#' @examples
#' # Example 1: Visualize fit for simulated outbreak data
#' set.seed(123)
#' # Simulate ICC intervals from mixed distribution
#' icc_data <- c(
#'   rnorm(20, mean = 0, sd = 2),      # Co-primary cases
#'   rnorm(50, mean = 12, sd = 3),     # Primary-secondary cases
#'   rnorm(20, mean = 24, sd = 4)      # Primary-tertiary cases
#' )
#' icc_data <- round(pmax(icc_data, 0))  # Ensure non-negative
#'
#' # Plot with estimated parameters
#' plot_si_fit(
#'   dat = icc_data,
#'   mean = 12.5,
#'   sd = 3.2,
#'   weights = c(0.2, 0.6, 0.15, 0.05),
#'   dist = "normal"
#' )
#'
#' # Example 2: Using gamma distribution
#' plot_si_fit(
#'   dat = icc_data,
#'   mean = 12.0,
#'   sd = 3.5,
#'   weights = c(0.25, 0.65, 0.10),
#'   dist = "gamma",
#'   scaling_factor = 0.8
#' )
#'
plot_si_fit <- function(
  dat,
  mean,
  sd,
  weights,
  dist = "normal",
  scaling_factor = 1
) {
  if (dist == "gamma") {
    # Define breaks
    breaks <- seq(min(dat) - 0.51, max(dat) + 0.51, by = 1)

    # Create the plot
    p <- ggplot(data = data.frame(x = dat), aes(x = .data$x)) +
      geom_histogram(
        aes(y = after_stat(.data$density) * scaling_factor),
        breaks = breaks,
        fill = "lightblue",
        color = "black"
      ) +
      stat_function(
        fun = f_gam,
        args = list(
          w1 = weights[1],
          w2 = weights[2],
          w3 = weights[3],
          mu = mean,
          sigma = sd
        ),
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Index-case to case interval (days)", y = "Probability") +
      theme_minimal()
  }

  if (dist == "normal") {
    # Define breaks for histogram
    breaks <- seq(min(dat) - 0.51, max(dat) + 0.51, by = 1)

    # Create the plot
    p <- ggplot(data = data.frame(x = dat), aes(x = .data$x)) +
      geom_histogram(
        aes(y = after_stat(.data$density) * scaling_factor),
        breaks = breaks,
        fill = "lightblue",
        color = "black"
      ) +
      stat_function(
        fun = f_norm,
        args = list(
          w1 = weights[1],
          w2 = weights[2],
          w3 = weights[3],
          mu = mean,
          sigma = sd
        ),
        color = "red",
        linetype = "solid",
        linewidth = 1
      ) +
      labs(x = "Index-case to case interval (days)", y = "Density") +
      theme_minimal() +
      geom_vline(xintercept = mean, linetype = "dashed", color = "black")
  }

  return(p)
}

#' Plot Serial Interval Fit from si_estim Result
#'
#' A convenience wrapper for \code{\link{plot_si_fit}} that accepts the output from
#' \code{\link{si_estim}} directly, automatically handling the weight aggregation
#' for different distribution types.
#'
#' @param si_result list; the output from \code{\link{si_estim}} containing mean, sd,
#'                  and wts (weights) components
#' @param dat numeric vector; the index case-to-case (ICC) intervals in days used
#'            for estimation
#' @param dist character; the distribution family used for estimation. Must be either
#'             "normal" (default) or "gamma". Should match the distribution used in
#'             the original \code{si_estim()} call
#' @param scaling_factor numeric; multiplicative factor to adjust the height of the
#'                       fitted density curve. Defaults to 1
#'
#' @return A \code{ggplot2} object showing the fitted distribution overlaid on a
#'         histogram of the observed data
#'
#' @details
#' This function simplifies the plotting workflow by automatically aggregating the

#' component weights from \code{si_estim()} output:
#' \itemize{
#'   \item For \strong{normal distribution}: Aggregates 7 weights into 4 transmission
#'         route weights (co-primary, primary-secondary, primary-tertiary,
#'         primary-quaternary)
#'   \item For \strong{gamma distribution}: Uses the 4 weights directly (components
#'         1, 2, 4, 6 mapped to the 4 transmission routes)
#' }
#'
#' @seealso \code{\link{si_estim}} for serial interval estimation,
#'          \code{\link{plot_si_fit}} for the underlying plotting function
#'
#' @export
#' @examples
#' # Simulate some ICC interval data
#' set.seed(123)
#' icc_data <- c(
#'   abs(rnorm(15, mean = 0, sd = 2)),
#'   rnorm(40, mean = 12, sd = 3),
#'   rnorm(15, mean = 24, sd = 4)
#' )
#' icc_data <- round(pmax(icc_data, 0))
#'
#' # Estimate serial interval
#' \donttest{
#' result <- si_estim(icc_data, n = 50)
#'
#' # Plot using the convenience wrapper
#' plot_si_fit_result(result, icc_data, dist = "normal")
#' }
#'
plot_si_fit_result <- function(
  si_result,
  dat,
  dist = c("normal", "gamma"),
  scaling_factor = 1
) {
  dist <- match.arg(dist)

 if (dist == "normal") {
    # Aggregate 7 weights into 4 for normal distribution
    # Component 1: Co-primary
    # Components 2+3: Primary-secondary
    # Components 4+5: Primary-tertiary
    # Components 6+7: Primary-quaternary
    weights <- c(
      si_result$wts[1],
      si_result$wts[2] + si_result$wts[3],
      si_result$wts[4] + si_result$wts[5],
      si_result$wts[6] + si_result$wts[7]
    )
  } else {
    # Gamma distribution uses 4 weights directly (components 1, 2, 4, 6)
    weights <- si_result$wts[1:3]
  }

  plot_si_fit(
    dat = dat,
    mean = si_result$mean,
    sd = si_result$sd,
    weights = weights,
    dist = dist,
    scaling_factor = scaling_factor
  )
}
