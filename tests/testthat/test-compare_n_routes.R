# Tests for compare_n_routes()

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
  local_mocked_bindings(
    si_estim = mock_si_estim,
    plot_si_fit_result = function(...) invisible(NULL)
  )
  capture.output(vis <- withVisible(compare_n_routes(1:20, n_routes_max = 3)))
  expect_false(vis$visible)
})

test_that("compare_n_routes rejects invalid n_routes_max", {
  for (bad in list(1, 0, -3, 2.5, NA, Inf, "4", c(3, 4))) {
    expect_error(compare_n_routes(1:20, n_routes_max = bad),
                 regexp = "integer >= 2")
  }
})

test_that("compare_n_routes accepts the smallest valid n_routes_max", {
  local_mocked_bindings(
    si_estim = mock_si_estim,
    plot_si_fit_result = function(...) invisible(NULL)
  )
  capture.output(res <- compare_n_routes(1:20, n_routes_max = 2))
  expect_equal(res$n_routes, 2L)
})

# -----------------------------------------------------------------------
# compare_n_routes: end-to-end with real fits
# -----------------------------------------------------------------------

test_that("compare_n_routes runs end to end on simulated data", {
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
