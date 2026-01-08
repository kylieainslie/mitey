# Create test file for si_estim function
test_that("si_estim produces correct estimates with simulated data", {
  library(fdrtool)
  set.seed(1234)

  # Simulate data with known parameters (as seen in your validation document)
  N <- 1000
  true_mu <- 15
  true_sigma <- 3
  hw1 <- 0.2
  hw2 <- 0.5
  hw3 <- 0.2
  hw4 <- 0.1

  # Generate from different transmission routes
  CP <- round(rhalfnorm((hw1*N), theta=sqrt(pi/2)/(sqrt(2)*true_sigma)))
  PS <- round(rnorm(hw2*N, mean=true_mu, sd=true_sigma))
  PT <- round(rnorm(hw3*N, mean=2*true_mu, sd=sqrt(2)*true_sigma))
  PQ <- round(rnorm(hw4*N, mean=3*true_mu, sd=sqrt(3)*true_sigma))

  sim_data <- c(CP, PS, PT, PQ)

  # Run si_estim with default parameters (normal distribution)
  result_normal <- si_estim(sim_data)

  # Test normal distribution estimates
  expect_type(result_normal, "list")
  expect_setequal(names(result_normal), c("mean", "sd", "wts", "converged", "iterations"))

  # Check that the estimate is close to true value
  expect_true(abs(result_normal$mean[1] - true_mu) < 1,
              "Mean estimate should be within 1 unit of true value")
  expect_true(abs(result_normal$sd[1] - true_sigma) < 1,
              "SD estimate should be within 1 unit of true value")

  # Check weights
  expect_equal(length(result_normal$wts), 7)  # Expected number of weight components
  expect_true(abs(result_normal$wts[1] - hw1) < 0.1,
              "Co-primary weight should be close to true value")
  expect_true(abs((result_normal$wts[2] + result_normal$wts[3]) - hw2) < 0.1,
              "Primary-secondary weight should be close to true value")
  expect_true(abs((result_normal$wts[4] + result_normal$wts[5]) - hw3) < 0.1,
              "Primary-tertiary weight should be close to true value")
  expect_true(abs((result_normal$wts[6] + result_normal$wts[7]) - hw4) < 0.1,
              "Primary-quaternary weight should be close to true value")

  # Test gamma distribution
  result_gamma <- si_estim(sim_data, dist = "gamma")

  # Basic validation for gamma distribution
  expect_type(result_gamma, "list")
  expect_named(result_gamma, c("mean", "sd", "wts", "converged", "iterations"))

  # Gamma may not match as closely since data was generated using normal distribution
  # But should still be reasonable
  expect_true(result_gamma$mean[1] > 5 && result_gamma$mean[1] < 25,
              "Gamma mean should be within reasonable range")
  expect_true(result_gamma$sd[1] > 0, "SD should be positive")
})

test_that("si_estim handles edge cases appropriately", {
  # Test with minimal dataset (at least 2 values with some difference)
  small_data <- c(1, 5, 10)
  result_small <- si_estim(small_data)
  expect_type(result_small, "list")
  expect_true(!is.na(result_small$mean[1]))
  expect_true(!is.na(result_small$sd[1]))

  # Test with initial values provided
  init_test <- si_estim(c(10, 12, 15, 20), init = c(12, 3))
  expect_type(init_test, "list")
  expect_length(init_test$mean, 1)
  expect_length(init_test$sd, 1)

  # Test with a dataset that has all the same values
  # This should still work but might give a warning about standard deviation
  repeated_data <- rep(10, 5)
  expect_error(
    si_estim(repeated_data),
    regexp = "Initial standard deviation must be positive|standard deviation|zero variance|singular",
    info = "Should error when all values are identical (SD = 0)"
  )

  if (exists("result_repeated")) {
    expect_type(result_repeated, "list")
    expect_equal(result_repeated$mean[1], 10)
    # SD might be 0 or close to 0
    expect_true(result_repeated$sd[1] >= 0)
  }
})

test_that("si_estim handles invalid inputs appropriately", {
  # Test with negative values (should now give a warning)
  expect_warning(
    result_negative <- si_estim(c(-5, 10, 15)),
    regexp = "negative values|negative serial intervals"
  )
  # Check that it still returns valid results despite the warning
  expect_type(result_negative, "list")
  expect_true(!is.na(result_negative$mean[1]))
  expect_true(!is.na(result_negative$sd[1]))

  # Test with NAs (should now throw an error)
  expect_error(
    si_estim(c(NA, 10, 15)),
    regexp = "Data contains NA values|NA values"
  )

  # Test with non-numeric data
  expect_error(
    si_estim(as.character(c(1, 2, 3))),
    regexp = "numeric|character|type"
  )

  # Test with invalid distribution
  expect_error(
    si_estim(c(10, 15, 20), dist = "invalid_dist"),
    regexp = "normal|gamma|invalid"
  )

  # Test with a single value (should error)
  expect_error(
    si_estim(10),
    regexp = "at least 2|single value|insufficient"
  )

  # Test gamma distribution with negative values (should error)
  expect_condition(
   expect_error(
      si_estim(c(-5, 10, 15), dist = "gamma"),
      regexp = "Gamma distribution cannot be used with negative values|negative values.*gamma"
    ),
    regexp = "negative values|negative serial intervals",
    class = "warning"
  )
})

test_that("si_estim properly estimates historical data", {
  # Create a small subset of data similar to what's in the validation document
  # This test is for simple functionality, not exact matching of historical results
  measles_data <- c(13, 14, 14, 14, 15, 15, 16, 16, 16, 18)
  influenza_data <- c(1, 2, 2, 3, 3, 3, 4, 4, 5)

  # Test with measles data
  measles_result <- si_estim(measles_data, init = c(15, 3))
  expect_type(measles_result, "list")
  expect_true(measles_result$mean[1] > 10 && measles_result$mean[1] < 20,
              "Measles mean should be in reasonable range")

  # Test with influenza data
  flu_result <- si_estim(influenza_data, init = c(3, 1))
  expect_type(flu_result, "list")
  expect_true(flu_result$mean[1] > 0 && flu_result$mean[1] < 10,
              "Influenza mean should be in reasonable range")
})

test_that("si_estim is consistent with different seeds", {
  set.seed(123)
  data1 <- rnorm(100, mean = 15, sd = 3)

  # Run with different seeds
  set.seed(456)
  result1 <- si_estim(data1)

  set.seed(789)
  result2 <- si_estim(data1)

  # Results should be the same since the data is the same
  expect_equal(result1$mean, result2$mean)
  expect_equal(result1$sd, result2$sd)
})

test_that("si_estim converges with different initial values", {
  set.seed(123)
  test_data <- rnorm(50, mean = 10, sd = 2)

  # Run with different initial values
  result1 <- si_estim(test_data, init = c(5, 1))
  result2 <- si_estim(test_data, init = c(15, 3))

  # Results should be similar regardless of initial values
  # We use a tolerance because EM algorithm might converge to slightly different values
  expect_equal(result1$mean, result2$mean, tolerance = 1)
  expect_equal(result1$sd, result2$sd, tolerance = 1)
})

test_that("si_estim convergence diagnostics work correctly", {
  set.seed(123)
  test_data <- rnorm(100, mean = 10, sd = 2)


  # Test that convergence info is returned
  result <- si_estim(test_data)
  expect_true("converged" %in% names(result))
  expect_true("iterations" %in% names(result))
  expect_type(result$converged, "logical")
  expect_type(result$iterations, "integer")

  # Test early stopping with default tolerance
  # Should converge before max iterations for well-behaved data
  result_default <- si_estim(test_data, n = 100)
  expect_true(result_default$iterations <= 100)


  # Test with tolerance = 0 (no early stopping - always runs all iterations)
  result_no_tol <- si_estim(test_data, n = 20, tol = 0)
  expect_equal(result_no_tol$iterations, 20)
  # When tol = 0, convergence checking is disabled, so converged is always FALSE
  expect_false(result_no_tol$converged)

  # Test with very loose tolerance (should converge quickly)
  result_loose <- si_estim(test_data, n = 100, tol = 0.1)
  expect_true(result_loose$converged)
  expect_true(result_loose$iterations < 100)

  # Test invalid tolerance
  expect_error(si_estim(test_data, tol = -1), "non-negative")
  expect_error(si_estim(test_data, tol = "invalid"), "numeric")
  expect_error(si_estim(test_data, tol = c(1, 2)), "single")
  expect_error(si_estim(test_data, tol = NA), "finite")
  expect_error(si_estim(test_data, tol = Inf), "finite")
  expect_error(si_estim(test_data, tol = NaN), "finite")
})
