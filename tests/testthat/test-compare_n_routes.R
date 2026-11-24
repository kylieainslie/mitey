# Tests for compare_n_routes(), n_optimal1() and n_optimal2()

# -----------------------------------------------------------------------
# n_optimal1 / n_optimal2
# -----------------------------------------------------------------------

test_that("n_optimal1 returns floor of 95th percentile / mean", {
  dat <- 1:10
  # quantile(1:10, 0.95) = 9.55, mean = 5.5 -> 1.736...
  expect_identical(n_optimal1(dat), 1L)

  dat2 <- c(rep(1, 19), 100)
  # quantile = 5.95, mean = 5.95 -> 1
  expect_identical(n_optimal1(dat2), 1L)

  dat3 <- c(rep(1, 10), rep(30, 10))
  # quantile = 30, mean = 15.5 -> 1.935...
  expect_identical(n_optimal1(dat3), 1L)

  dat4 <- c(rep(0, 18), 50, 50)
  # quantile = 50, mean = 5 -> 10
  expect_identical(n_optimal1(dat4), 10L)
})

test_that("n_optimal2 returns floor of max / sd", {
  dat <- 1:10
  # 10 / sd(1:10) = 10 / 3.0277 = 3.30...
  expect_identical(n_optimal2(dat), 3L)

  dat2 <- c(0, 10)
  # 10 / sd(c(0, 10)) = 10 / 7.071 = 1.41...
  expect_identical(n_optimal2(dat2), 1L)

  set.seed(1)
  dat3 <- rnorm(200, mean = 20, sd = 2)
  expect_identical(n_optimal2(dat3), as.integer(max(dat3) / sd(dat3)))
})

test_that("n_optimal heuristics return a single integer", {
  set.seed(42)
  dat <- round(abs(c(rnorm(30, 0, 2), rnorm(80, 12, 3), rnorm(30, 24, 4))))
  for (res in list(n_optimal1(dat), n_optimal2(dat))) {
    expect_type(res, "integer")
    expect_length(res, 1)
  }
})

# -----------------------------------------------------------------------
# compare_n_routes: logic checked against a mocked si_estim()
# -----------------------------------------------------------------------

# Fake fit with fixed log-likelihoods per n_routes, so AIC/BIC are known.
fake_loglik <- c(`2` = -100, `3` = -90, `4` = -89)

mock_si_estim <- function(dat, n = 50, dist = "normal", n_routes = 4L, ...) {
  list(
    mean       = 10 + n_routes,
    sd         = 2,
    wts        = rep(1 / n_routes, n_routes - 1),
    converged  = TRUE,
    iterations = 7L,
    loglik     = unname(fake_loglik[as.character(n_routes)])
  )
}

test_that("compare_n_routes computes parameters, AIC and BIC (normal)", {
  skip_if_not_installed("dplyr")
  local_mocked_bindings(
    si_estim = mock_si_estim,
    plot_si_fit_result = function(...) invisible(NULL)
  )
  dat <- 1:20

  expect_output(
    res <- compare_n_routes(dat, n_routes_max = 4, dist = "normal"),
    "Comparison of the n_routes"
  )

  expect_s3_class(res, "data.frame")
  expect_named(res, c("n_routes", "mu", "sigma", "loglik", "aic", "bic",
                      "n_params", "converged", "iterations"))
  expect_equal(res$n_routes, 2:4)

  # Normal: p = 2 + (2 * n_routes - 2)
  expect_equal(res$n_params, c(4L, 6L, 8L))
  expect_equal(res$loglik, unname(fake_loglik))
  expect_equal(res$aic, -2 * res$loglik + 2 * res$n_params)
  expect_equal(res$bic, -2 * res$loglik + log(length(dat)) * res$n_params)
  expect_equal(res$mu, 12:14)
  expect_true(all(res$converged))
  expect_equal(res$iterations, rep(7L, 3))
})

test_that("compare_n_routes uses gamma parameter count", {
  skip_if_not_installed("dplyr")
  local_mocked_bindings(
    si_estim = mock_si_estim,
    plot_si_fit_result = function(...) invisible(NULL)
  )

  expect_output(
    res <- compare_n_routes(1:20, n_routes_max = 4, dist = "gamma")
  )

  # Gamma: p = 2 + (n_routes - 1)
  expect_equal(res$n_params, c(3L, 4L, 5L))
  expect_equal(res$aic, -2 * res$loglik + 2 * res$n_params)
})

test_that("compare_n_routes reports the model with the lowest AIC and BIC", {
  skip_if_not_installed("dplyr")
  local_mocked_bindings(
    si_estim = mock_si_estim,
    plot_si_fit_result = function(...) invisible(NULL)
  )

  # n = 20, normal: AIC = 208, 192, 194; BIC = 212.0, 198.0, 202.0 -> both 3
  out <- capture.output(res <- compare_n_routes(1:20, n_routes_max = 4))
  expect_true(any(grepl("according to BIC : 3", out)))
  expect_true(any(grepl("according to AIC : 3", out)))
  expect_equal(res$n_routes[which.min(res$bic)], 3L)
})

test_that("compare_n_routes forwards arguments to si_estim", {
  skip_if_not_installed("dplyr")
  calls <- list()
  local_mocked_bindings(
    si_estim = function(dat, n, dist, n_routes, ...) {
      calls[[length(calls) + 1]] <<- list(n = n, dist = dist,
                                          n_routes = n_routes, extra = list(...))
      mock_si_estim(dat, n = n, dist = dist, n_routes = n_routes)
    },
    plot_si_fit_result = function(...) invisible(NULL)
  )

  expect_output(
    compare_n_routes(1:20, n_routes_max = 3, dist = "gamma", n = 12,
                     wind = 2, n_starts = 3)
  )

  expect_length(calls, 2)
  expect_equal(vapply(calls, `[[`, integer(1), "n_routes"), 2:3)
  for (cl in calls) {
    expect_equal(cl$n, 12)
    expect_equal(cl$dist, "gamma")
    expect_equal(cl$extra, list(wind = 2, n_starts = 3))
  }
})

test_that("compare_n_routes returns its results invisibly", {
  skip_if_not_installed("dplyr")
  local_mocked_bindings(
    si_estim = mock_si_estim,
    plot_si_fit_result = function(...) invisible(NULL)
  )
  capture.output(vis <- withVisible(compare_n_routes(1:20, n_routes_max = 3)))
  expect_false(vis$visible)
})

# -----------------------------------------------------------------------
# compare_n_routes: end-to-end with real fits
# -----------------------------------------------------------------------

test_that("compare_n_routes runs end to end on simulated data", {
  skip_if_not_installed("dplyr")
  skip_on_cran()
  set.seed(42)
  icc <- round(pmax(c(abs(rnorm(30, 0, 2)), rnorm(80, 12, 3),
                      rnorm(30, 24, 4)), 0))

  # Plots are printed as a side effect; render them so that errors inside
  # ggplot layers (which ggplot reports only as warnings) are caught.
  # Other warnings, e.g. the gamma density being undefined at x = 0, are fine.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (dist in c("normal", "gamma")) {
    warns <- character()
    withCallingHandlers(
      expect_output(
        res <- compare_n_routes(icc, n_routes_max = 3, dist = dist, n = 20),
        "Best n_routes according to BIC"
      ),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_false(any(grepl("Computation failed", warns)), label = dist)
    expect_equal(nrow(res), 2)
    expect_equal(res$n_routes, 2:3)
    expect_true(all(is.finite(res$loglik)))
    expect_true(all(is.finite(res$aic)))
    expect_true(all(is.finite(res$bic)))
    expect_true(all(res$mu > 0))
    expect_true(all(res$sigma > 0))
    expect_type(res$converged, "logical")
  }
})
