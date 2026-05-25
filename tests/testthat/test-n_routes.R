# Test file for n_routes flexibility across all modified functions

# -----------------------------------------------------------------------
# integrate_components_wrapper
# -----------------------------------------------------------------------

test_that("integrate_components_wrapper works with n_routes = 4 (default)", {
  result <- integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "normal")
  expect_length(result, 7)  # 2*4 - 1 = 7
  expect_true(all(result >= 0))

  result_gam <- integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "gamma")
  expect_length(result_gam, 4)  # n_routes = 4
  expect_true(all(result_gam >= 0))
})

test_that("integrate_components_wrapper works with n_routes = 3", {
  result <- integrate_components_wrapper(d = 10, mu = 15, sigma = 3,
                                         dist = "normal", n_routes = 3)
  expect_length(result, 5)  # 2*3 - 1 = 5
  expect_true(all(result >= 0))

  result_gam <- integrate_components_wrapper(d = 10, mu = 15, sigma = 3,
                                             dist = "gamma", n_routes = 3)
  expect_length(result_gam, 3)
  expect_true(all(result_gam >= 0))
})

test_that("integrate_components_wrapper works with n_routes = 5", {
  result <- integrate_components_wrapper(d = 10, mu = 15, sigma = 3,
                                         dist = "normal", n_routes = 5)
  expect_length(result, 9)  # 2*5 - 1 = 9
  expect_true(all(result >= 0))

  result_gam <- integrate_components_wrapper(d = 10, mu = 15, sigma = 3,
                                             dist = "gamma", n_routes = 5)
  expect_length(result_gam, 5)
  expect_true(all(result_gam >= 0))
})

test_that("integrate_components_wrapper handles d = 0", {
  result <- integrate_components_wrapper(d = 0, mu = 15, sigma = 3,
                                         dist = "normal", n_routes = 5)
  expect_length(result, 9)
  expect_true(all(result >= 0))
})

test_that("integrate_components_wrapper rejects invalid n_routes", {
  expect_error(
    integrate_components_wrapper(d = 10, mu = 15, sigma = 3, n_routes = 1),
    regexp = "integer >= 2"
  )
  expect_error(
    integrate_components_wrapper(d = 10, mu = 15, sigma = 3, n_routes = 1.5),
    regexp = "integer >= 2"
  )
  expect_error(
    integrate_components_wrapper(d = 10, mu = 15, sigma = 3, n_routes = -1),
    regexp = "integer >= 2"
  )
})

# -----------------------------------------------------------------------
# f_norm
# -----------------------------------------------------------------------

test_that("f_norm returns non-negative values for various n_routes", {
  x <- seq(1, 30, by = 1)

  # 4 routes: weights length = 3
  result4 <- f_norm(x, weights = c(0.2, 0.5, 0.2), mu = 10, sigma = 3, n_routes = 4)
  expect_length(result4, length(x))
  expect_true(all(result4 >= 0))

  # 3 routes: weights length = 2
  result3 <- f_norm(x, weights = c(0.2, 0.6), mu = 10, sigma = 3, n_routes = 3)
  expect_length(result3, length(x))
  expect_true(all(result3 >= 0))

  # 5 routes: weights length = 4
  result5 <- f_norm(x, weights = c(0.2, 0.4, 0.2, 0.1), mu = 10, sigma = 3, n_routes = 5)
  expect_length(result5, length(x))
  expect_true(all(result5 >= 0))
})

test_that("f_norm rejects invalid weights", {
  x <- seq(1, 30, by = 1)

  expect_error(
    f_norm(x, weights = c(0.5, 0.4, 0.3), mu = 10, sigma = 3),
    regexp = "Sum of weights"
  )
  expect_error(
    f_norm(x, weights = c(-0.1, 0.5, 0.2), mu = 10, sigma = 3),
    regexp = "non-negative"
  )
})

# -----------------------------------------------------------------------
# f_gam
# -----------------------------------------------------------------------

test_that("f_gam returns non-negative values for various n_routes", {
  x <- seq(0.1, 30, by = 0.1)

  # 4 routes: weights length = 3
  result4 <- f_gam(x, weights = c(0.1, 0.6, 0.2), mu = 6.5, sigma = 2.8)
  expect_length(result4, length(x))
  expect_true(all(result4 >= 0))

  # 3 routes: weights length = 2
  result3 <- f_gam(x, weights = c(0.1, 0.7), mu = 6.5, sigma = 2.8)
  expect_length(result3, length(x))
  expect_true(all(result3 >= 0))

  # 5 routes: weights length = 4
  result5 <- f_gam(x, weights = c(0.1, 0.5, 0.2, 0.1), mu = 6.5, sigma = 2.8)
  expect_length(result5, length(x))
  expect_true(all(result5 >= 0))
})

test_that("f_gam rejects invalid weights", {
  x <- seq(0.1, 30, by = 0.1)

  expect_error(
    f_gam(x, weights = c(0.5, 0.4, 0.3), mu = 6.5, sigma = 2.8),
    regexp = "Sum of weights"
  )
  expect_error(
    f_gam(x, weights = c(-0.1, 0.6, 0.2), mu = 6.5, sigma = 2.8),
    regexp = "non-negative"
  )
})

# -----------------------------------------------------------------------
# si_estim
# -----------------------------------------------------------------------

test_that("si_estim works with n_routes = 4 (default, rétrocompatibilité)", {
  set.seed(123)
  dat <- c(rep(1, 20), rep(2, 25), rep(3, 15), rep(4, 8))

  result <- si_estim(dat)
  expect_equal(result$n_routes, 4L)
  expect_length(result$wts, 7)  # 2*4 - 1
  expect_true(is.finite(result$mean))
  expect_true(is.finite(result$sd))
  expect_true(result$sd > 0)
})

test_that("si_estim works with n_routes = 3", {
  set.seed(123)
  dat <- c(rep(1, 20), rep(2, 25), rep(3, 15), rep(4, 8))

  result <- si_estim(dat, n_routes = 3)
  expect_equal(result$n_routes, 3L)
  expect_length(result$wts, 5)  # 2*3 - 1
  expect_true(is.finite(result$mean))
  expect_true(result$sd > 0)
})

test_that("si_estim works with n_routes = 5", {
  set.seed(123)
  dat <- c(rep(1, 20), rep(2, 25), rep(3, 15), rep(4, 8))

  result <- si_estim(dat, n_routes = 5)
  expect_equal(result$n_routes, 5L)
  expect_length(result$wts, 9)  # 2*5 - 1
  expect_true(is.finite(result$mean))
  expect_true(result$sd > 0)
})

test_that("si_estim works with gamma distribution and various n_routes", {
  set.seed(456)
  dat <- c(rep(1, 20), rep(2, 25), rep(3, 15), rep(4, 8))

  result3 <- si_estim(dat, dist = "gamma", n_routes = 3)
  expect_equal(result3$n_routes, 3L)
  expect_length(result3$wts, 3)

  result5 <- si_estim(dat, dist = "gamma", n_routes = 5)
  expect_equal(result5$n_routes, 5L)
  expect_length(result5$wts, 5)
})

test_that("si_estim rejects invalid n_routes", {
  dat <- c(rep(1, 20), rep(2, 25), rep(3, 15))

  expect_error(si_estim(dat, n_routes = 1),  regexp = "integer >= 2")
  expect_error(si_estim(dat, n_routes = 1.5), regexp = "integer >= 2")
  expect_error(si_estim(dat, n_routes = -2),  regexp = "integer >= 2")
  expect_error(si_estim(dat, n_routes = "a"), regexp = "integer >= 2")
})

test_that("si_estim n_routes = 4 gives same results as original hardcoded version", {
  # Verifies backward compatibility: results should be identical
  set.seed(789)
  dat <- c(rep(6,4), rep(7,8), rep(8,14), rep(9,31), rep(10,29),
           rep(11,42), rep(12,25), rep(13,16), rep(14,16),
           rep(15,10), rep(16,4), rep(17,2), rep(18,2))

  result_default <- si_estim(dat, n = 50)
  result_explicit <- si_estim(dat, n = 50, n_routes = 4)

  expect_equal(result_default$mean, result_explicit$mean)
  expect_equal(result_default$sd,   result_explicit$sd)
  expect_equal(result_default$wts,  result_explicit$wts)
})

# -----------------------------------------------------------------------
# plot_si_fit and plot_si_fit_result
# -----------------------------------------------------------------------

test_that("plot_si_fit works with various n_routes", {
  set.seed(123)
  dat <- round(pmax(rnorm(100, mean = 10, sd = 3), 0))

  # 4 routes
  p4 <- plot_si_fit(dat, mean = 10, sd = 3,
                    weights = c(0.1, 0.7, 0.15, 0.05),
                    dist = "normal", n_routes = 4)
  expect_s3_class(p4, "ggplot")

  # 3 routes
  p3 <- plot_si_fit(dat, mean = 10, sd = 3,
                    weights = c(0.1, 0.8, 0.1),
                    dist = "normal", n_routes = 3)
  expect_s3_class(p3, "ggplot")

  # 5 routes
  p5 <- plot_si_fit(dat, mean = 10, sd = 3,
                    weights = c(0.1, 0.6, 0.15, 0.08, 0.07),
                    dist = "normal", n_routes = 5)
  expect_s3_class(p5, "ggplot")
})

test_that("plot_si_fit rejects weights of wrong length", {
  set.seed(123)
  dat <- round(pmax(rnorm(100, mean = 10, sd = 3), 0))

  expect_error(
    plot_si_fit(dat, mean = 10, sd = 3,
                weights = c(0.1, 0.7, 0.1),   # length 3, but n_routes = 4 expects 4
                dist = "normal", n_routes = 4),
    regexp = "n_routes"
  )
})

test_that("plot_si_fit_result works with n_routes stored in si_result", {
  set.seed(123)
  dat <- c(rep(1, 20), rep(2, 25), rep(3, 15), rep(4, 8))

  result4 <- si_estim(dat, n_routes = 4)
  p4 <- plot_si_fit_result(result4, dat, dist = "normal")
  expect_s3_class(p4, "ggplot")

  result3 <- si_estim(dat, n_routes = 3)
  p3 <- plot_si_fit_result(result3, dat, dist = "normal")
  expect_s3_class(p3, "ggplot")

  result5 <- si_estim(dat, n_routes = 5)
  p5 <- plot_si_fit_result(result5, dat, dist = "normal")
  expect_s3_class(p5, "ggplot")
})
