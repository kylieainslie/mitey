#' Visualize Serial Interval Distribution Fit to Outbreak Data
#'
#' Creates a diagnostic plot showing the fitted serial interval mixture distribution
#' overlaid on a histogram of observed index case-to-case (ICC) intervals from outbreak
#' data.
#'
#' @param dat numeric vector; the index case-to-case (ICC) intervals in days.
#' @param mean numeric; the estimated mean of the serial interval distribution in days.
#' @param sd numeric; the estimated standard deviation of the serial interval
#'           distribution in days.
#' @param weights numeric vector of length n_routes - 1; the estimated weights for
#'   each transmission route component, starting from Co-primary. The last route
#'   weight is derived internally as 1 - sum(weights).
#' @param dist character; the distribution family. Must be either "normal" (default)
#'   or "gamma".
#' @param scaling_factor numeric; multiplicative factor to adjust the height of the
#'   fitted density curve. Defaults to 1.
#' @param n_routes integer; number of transmission routes modelled. Must be >= 2.
#'   Defaults to 4. Must match the value used in \code{si_estim()}.
#'
#' @return A \code{ggplot2} object.
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
#' # Example 1: 4 routes, normal distribution
#' set.seed(123)
#' icc_data <- c(
#'   rnorm(20, mean = 0, sd = 2),
#'   rnorm(50, mean = 12, sd = 3),
#'   rnorm(20, mean = 24, sd = 4)
#' )
#' icc_data <- round(pmax(icc_data, 0))
#'
#' plot_si_fit(
#'   dat = icc_data,
#'   mean = 12.5,
#'   sd = 3.2,
#'   weights = c(0.2, 0.6, 0.15),
#'   dist = "normal",
#'   n_routes = 4
#' )
#'
#' # Example 2: 5 routes, gamma distribution
#' plot_si_fit(
#'   dat = icc_data,
#'   mean = 12.0,
#'   sd = 3.5,
#'   weights = c(0.25, 0.65, 0.05, 0.03),
#'   dist = "gamma",
#'   n_routes = 5,
#'   scaling_factor = 0.8
#' )
#'
plot_si_fit <- function(
    dat,
    mean,
    sd,
    weights,
    dist           = "normal",
    scaling_factor = 1,
    n_routes       = 4L
) {
  n_routes <- as.integer(n_routes)

  if (dist == "normal") {
    if (length(weights) != 2L * n_routes - 1L) {
      stop(paste0(
        "weights must have length 2*n_routes - 1 = ", 2L * n_routes - 1L,
        ", but has length ", length(weights), "."
      ))
    }
  } else {
    if (length(weights) != n_routes - 1L) {
      stop(paste0(
        "weights must have length n_routes - 1 = ", n_routes - 1L,
        ", but has length ", length(weights), "."
      ))
    }
  }

  breaks <- seq(min(dat) - 0.51, max(dat) + 0.51, by = 1)

  if (dist == "gamma") {
    p <- ggplot(data = data.frame(x = dat), aes(x = .data$x)) +
      geom_histogram(
        aes(y = after_stat(.data$density) * scaling_factor),
        breaks = breaks,
        fill   = "lightblue",
        color  = "black"
      ) +
      stat_function(
        fun  = f_gam,
        args = list(
          weights = weights,
          mu      = mean,
          sigma   = sd,
          n_routes = n_routes
        ),
        color     = "red",
        linetype  = "dashed",
        linewidth = 1
      ) +
      labs(x = "Index-case to case interval (days)", y = "Probability") +
      theme_minimal()
  }

  if (dist == "normal") {
    p <- ggplot(data = data.frame(x = dat), aes(x = .data$x)) +
      geom_histogram(
        aes(y = after_stat(.data$density) * scaling_factor),
        breaks = breaks,
        fill   = "lightblue",
        color  = "black"
      ) +
      stat_function(
        fun  = f_norm,
        args = list(
          weights  = weights,
          mu       = mean,
          sigma    = sd,
          n_routes = n_routes
        ),
        color     = "red",
        linetype  = "solid",
        linewidth = 1
      ) +
      labs(
        x     = "Index-case to case interval (days)",
        y     = "Density",
        title = paste("Serial interval mixture fit -", n_routes, "transmission routes")
      ) +
      theme_minimal() +
      geom_vline(xintercept = mean, linetype = "dashed", color = "black")
  }

  return(p)
}


#' Plot Serial Interval Fit from si_estim Result
#'
#' A convenience wrapper for \code{\link{plot_si_fit}} that accepts the output from
#' \code{\link{si_estim}} directly, automatically handling the weight aggregation
#' for different distribution types and number of routes.
#'
#' @param si_result list; the output from \code{\link{si_estim}} containing mean, sd,
#'   wts, and n_routes components.
#' @param dat numeric vector; the index case-to-case (ICC) intervals in days.
#' @param dist character; the distribution family. Must be either "normal" (default)
#'   or "gamma". Should match the distribution used in \code{si_estim()}.
#' @param scaling_factor numeric; multiplicative factor to adjust the height of the
#'   fitted density curve. Defaults to 1.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' This function reads \code{n_routes} directly from the \code{si_result} object and
#' aggregates component weights automatically:
#' \itemize{
#'   \item For normal distribution: the Co-primary weight is taken directly, and
#'     weights for routes 2 to n_routes are aggregated by summing the two symmetric
#'     components (positive and negative).
#'   \item For gamma distribution: weights are taken directly as the first
#'     n_routes - 1 values from \code{si_result$wts}.
#' }
#'
#' @seealso \code{\link{si_estim}}, \code{\link{plot_si_fit}}
#'
#' @export
#' @examples
#' set.seed(123)
#' icc_data <- c(
#'   abs(rnorm(15, mean = 0, sd = 2)),
#'   rnorm(40, mean = 12, sd = 3),
#'   rnorm(15, mean = 24, sd = 4)
#' )
#' icc_data <- round(pmax(icc_data, 0))
#'
#' \donttest{
#' # 4 routes (default)
#' result <- si_estim(icc_data, n = 50)
#' plot_si_fit_result(result, icc_data, dist = "normal")
#'
#' # 5 routes
#' result5 <- si_estim(icc_data, n = 50, n_routes = 5)
#' plot_si_fit_result(result5, icc_data, dist = "normal")
#' }
#'
plot_si_fit_result <- function(
    si_result,
    dat,
    dist           = c("normal", "gamma"),
    scaling_factor = 1
) {
  dist     <- match.arg(dist)
  n_routes <- as.integer(si_result$n_routes)

  if (dist == "normal") {
    weights <- si_result$wts  # poids bruts, longueur 2*n_routes - 1
  } else {
    weights <- si_result$wts[seq_len(n_routes - 1L)]
  }


  plot_si_fit(
    dat            = dat,
    mean           = si_result$mean,
    sd             = si_result$sd,
    weights        = weights,
    dist           = dist,
    scaling_factor = scaling_factor,
    n_routes       = n_routes
  )
}
