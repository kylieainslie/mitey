
# Test helpers
create_test_data <- function() {
  # Create synthetic data for testing
  dates <- as.Date("2023-01-01") + 0:9
  incidence <- c(1, 2, 4, 8, 12, 15, 10, 6, 3, 1)
  list(
    dates = dates,
    incidence = incidence,
    si_mean = 3,
    si_sd = 1.5,
    n = length(dates)
  )
}

# Test the create_day_diff_matrix function
test_that("create_day_diff_matrix calculates correct day differences", {
  test_data <- create_test_data()
  day_diffs <- create_day_diff_matrix(test_data$dates)

  # Check dimensions
  expect_equal(dim(day_diffs), c(test_data$n, test_data$n))

  # Check diagonal is zero (same day)
  expect_equal(diag(day_diffs), rep(0, test_data$n))

  # Check some specific values
  expect_equal(day_diffs[2, 1], 1)  # Second date minus first date = 1 day
  expect_equal(day_diffs[1, 2], -1) # First date minus second date = -1 day
  expect_equal(day_diffs[5, 1], 4)  # Fifth date minus first date = 4 days
})

# Test the calculate_si_probability_matrix function
test_that("calculate_si_probability_matrix creates valid probability matrix", {
  test_data <- create_test_data()
  day_diffs <- create_day_diff_matrix(test_data$dates)

  # Test with gamma distribution
  si_prob_gamma <- calculate_si_probability_matrix(
    day_diffs, test_data$si_mean, test_data$si_sd, "gamma"
  )

  # Check dimensions
  expect_equal(dim(si_prob_gamma), c(test_data$n, test_data$n))

  # Check all values are non-negative
  expect_true(all(si_prob_gamma >= 0))

  # Check backward transmission is zero
  expect_equal(sum(si_prob_gamma[upper.tri(si_prob_gamma)]), 0)

  # Test with normal distribution
  si_prob_normal <- calculate_si_probability_matrix(
    day_diffs, test_data$si_mean, test_data$si_sd, "normal"
  )

  # Check dimensions
  expect_equal(dim(si_prob_normal), c(test_data$n, test_data$n))

  # Check all values are non-negative
  expect_true(all(si_prob_normal >= 0))

  # Check backward transmission is zero
  expect_equal(sum(si_prob_normal[upper.tri(si_prob_normal)]), 0)

  # Verify that values roughly match expected densities
  # For a positive time difference = 3 (which should match the mean)
  expected_idx <- which(day_diffs == 3, arr.ind = TRUE)[1, ]

  # For gamma distribution
  shape <- (test_data$si_mean / test_data$si_sd)^2
  rate <- test_data$si_mean / (test_data$si_sd^2)
  expected_density_gamma <- dgamma(3, shape = shape, rate = rate)

  # For normal distribution
  expected_density_normal <- dnorm(3, mean = test_data$si_mean, sd = test_data$si_sd)

  # Allow some numerical precision difference
  expect_equal(si_prob_gamma[expected_idx[1], expected_idx[2]], expected_density_gamma, tolerance = 1e-7)
  expect_equal(si_prob_normal[expected_idx[1], expected_idx[2]], expected_density_normal, tolerance = 1e-7)
})

# Test the smooth_estimates function
test_that("smooth_estimates correctly applies moving average smoothing", {
  # Create a test vector
  test_values <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

  # Test with window = 3
  smoothed_3 <- smooth_estimates(test_values, 3)
  expect_equal(smoothed_3[2], mean(c(1, 2, 3)))  # First complete window
  expect_equal(smoothed_3[5], mean(c(4, 5, 6)))  # Middle window
  expect_equal(smoothed_3[9], mean(c(8, 9, 10))) # Last complete window

  # Test with window = 5
  smoothed_5 <- smooth_estimates(test_values, 5)
  expect_equal(smoothed_5[3], mean(c(1, 2, 3, 4, 5)))  # First complete window
  expect_equal(smoothed_5[8], mean(c(6, 7, 8, 9, 10))) # Last complete window

  # Test with NA values
  test_values_na <- c(1, NA, 3, 4, NA, 6)
  smoothed_na <- smooth_estimates(test_values_na, 3)
  expect_equal(smoothed_na[3], mean(c(3, 4), na.rm = TRUE))

  # Test with Inf values
  test_values_inf <- c(1, 2, Inf, 4, 5)
  smoothed_inf <- smooth_estimates(test_values_inf, 3)
  expect_equal(smoothed_inf[1], mean(c(1, 2)))  # Inf should be excluded
  expect_equal(smoothed_inf[5], mean(c(4, 5)))  # Inf should be excluded
})

# Test the calculate_truncation_correction function
test_that("calculate_truncation_correction produces valid correction factors", {
  test_data <- create_test_data()

  # Test with gamma distribution
  correction_gamma <- calculate_truncation_correction(
    test_data$dates, test_data$si_mean, test_data$si_sd, "gamma"
  )

  # Check length
  expect_equal(length(correction_gamma), test_data$n)

  # Check early dates (should have corrections close to 1)
  expect_true(all(correction_gamma[1:5] >= 1 | is.na(correction_gamma[1:5])))

  # Test with normal distribution
  correction_normal <- calculate_truncation_correction(
    test_data$dates, test_data$si_mean, test_data$si_sd, "normal"
  )

  # Check length
  expect_equal(length(correction_normal), test_data$n)

  # Check early dates (should have corrections close to 1)
  expect_true(all(correction_normal[1:5] >= 1 | is.na(correction_normal[1:5])))
})

# Test the generate_case_bootstrap function
test_that("generate_case_bootstrap produces valid bootstrap samples", {
  # Create test incidence data
  incidence <- c(5, 10, 15, 20, 10)
  total_cases <- sum(incidence)

  # Generate bootstrap sample
  set.seed(123) # For reproducibility
  bootstrap_sample <- generate_case_bootstrap(incidence)

  # Check length
  expect_equal(length(bootstrap_sample), length(incidence))

  # Check total cases preserved
  expect_equal(sum(bootstrap_sample), total_cases)

  # Check all values are non-negative integers
  expect_true(all(bootstrap_sample >= 0))
  expect_true(all(bootstrap_sample == floor(bootstrap_sample)))

  # Check that running it twice gives different results (stochastic)
  set.seed(456) # Different seed
  bootstrap_sample2 <- generate_case_bootstrap(incidence)
  expect_false(identical(bootstrap_sample, bootstrap_sample2))
})

# Test the calculate_r_estimates function
test_that("calculate_r_estimates produces valid R estimates", {
  test_data <- create_test_data()
  day_diffs <- create_day_diff_matrix(test_data$dates)
  si_prob <- calculate_si_probability_matrix(
    day_diffs, test_data$si_mean, test_data$si_sd, "gamma"
  )

  # Calculate R without smoothing
  r_estimates <- calculate_r_estimates(
    test_data$incidence, si_prob, test_data$dates,
    test_data$si_mean, test_data$si_sd, "gamma", 0
  )

  # Check structure
  expect_true("r" %in% names(r_estimates))
  expect_true("r_corrected" %in% names(r_estimates))

  # Check length
  expect_equal(length(r_estimates$r), test_data$n)
  expect_equal(length(r_estimates$r_corrected), test_data$n)

  # For exponential growth with R = 2, early R estimates should be around 2
  # This is approximate since we're using synthetic data
  early_r_est <- mean(r_estimates$r[2:4], na.rm = TRUE)
  expect_true(early_r_est > 1 && early_r_est < 3.1)

  # Calculate R with smoothing
  r_estimates_smooth <- calculate_r_estimates(
    test_data$incidence, si_prob, test_data$dates,
    test_data$si_mean, test_data$si_sd, "gamma", 3
  )

  # Check smoothing was applied (values should be different)
  expect_false(identical(r_estimates$r, r_estimates_smooth$r))
})

# Test the calculate_bootstrap_ci function
test_that("calculate_bootstrap_ci produces valid confidence intervals", {
  test_data <- create_test_data()
  day_diffs <- create_day_diff_matrix(test_data$dates)
  si_prob <- calculate_si_probability_matrix(
    day_diffs, test_data$si_mean, test_data$si_sd, "gamma"
  )

  # Set smaller bootstrap samples for testing
  n_bootstrap <- 10

  # Calculate bootstrap CIs
  set.seed(123) # For reproducibility
  ci_estimates <- calculate_bootstrap_ci(
    test_data$incidence, si_prob, test_data$dates,
    test_data$si_mean, test_data$si_sd, "gamma", 0,
    n_bootstrap, 0.95
  )

  # Check structure
  expect_true(all(c("r_lower", "r_upper", "r_corrected_lower", "r_corrected_upper") %in% names(ci_estimates)))

  # Check length
  expect_equal(length(ci_estimates$r_lower), test_data$n)
  expect_equal(length(ci_estimates$r_upper), test_data$n)

  # Check lower bounds are less than or equal to upper bounds
  expect_true(all(ci_estimates$r_lower <= ci_estimates$r_upper |
                    (is.na(ci_estimates$r_lower) & is.na(ci_estimates$r_upper))))
  expect_true(all(ci_estimates$r_corrected_lower <= ci_estimates$r_corrected_upper |
                    (is.na(ci_estimates$r_corrected_lower) & is.na(ci_estimates$r_corrected_upper))))
})

# Test the main wallinga_lipsitch function
test_that("wallinga_lipsitch produces correct output format", {
  test_data <- create_test_data()

  # Basic run without bootstrap or shift
  result <- wallinga_lipsitch(
    test_data$incidence, test_data$dates,
    test_data$si_mean, test_data$si_sd
  )

  # Check columns
  expect_true(all(c("date", "incidence", "R", "R_corrected") %in% names(result)))
  expect_false("shifted_date" %in% names(result))
  expect_false("R_lower" %in% names(result))

  # Check dimensions
  expect_equal(nrow(result), test_data$n)

  # Run with bootstrap
  result_boot <- wallinga_lipsitch(
    test_data$incidence, test_data$dates,
    test_data$si_mean, test_data$si_sd,
    bootstrap = TRUE, n_bootstrap = 10
  )

  # Check additional columns
  expect_true(all(c("R_lower", "R_upper", "R_corrected_lower", "R_corrected_upper") %in% names(result_boot)))

  # Run with shift
  result_shift <- wallinga_lipsitch(
    test_data$incidence, test_data$dates,
    test_data$si_mean, test_data$si_sd,
    shift = TRUE
  )

  # Check shifted_date column
  expect_true("shifted_date" %in% names(result_shift))
  expect_equal(result_shift$shifted_date[1], result_shift$date[1] + round(test_data$si_mean))

  # Run with both bootstrap and shift
  result_both <- wallinga_lipsitch(
    test_data$incidence, test_data$dates,
    test_data$si_mean, test_data$si_sd,
    bootstrap = TRUE, n_bootstrap = 10,
    shift = TRUE
  )

  # Check all additional columns
  expect_true(all(c("R_lower", "R_upper", "R_corrected_lower", "R_corrected_upper", "shifted_date") %in% names(result_both)))
})

# Test edge cases and error handling
test_that("functions handle edge cases correctly", {
  test_data <- create_test_data()

  # Empty data
  expect_error(
    wallinga_lipsitch(numeric(0), as.Date(character(0)), test_data$si_mean, test_data$si_sd)
  )

  # Negative incidence
  expect_error(
    wallinga_lipsitch(c(-1, 1, 2), test_data$dates[1:3], test_data$si_mean, test_data$si_sd)
  )

  # Negative serial interval parameters
  expect_error(
    wallinga_lipsitch(test_data$incidence, test_data$dates, -1, test_data$si_sd)
  )
  expect_error(
    wallinga_lipsitch(test_data$incidence, test_data$dates, test_data$si_mean, -1)
  )

  # Invalid distribution
  expect_error(
    wallinga_lipsitch(test_data$incidence, test_data$dates, test_data$si_mean, test_data$si_sd, si_dist = "invalid")
  )

  # Data with all zeros
  zero_incidence <- rep(0, length(test_data$dates))
  result_zeros <- wallinga_lipsitch(zero_incidence, test_data$dates, test_data$si_mean, test_data$si_sd)
  expect_true(all(is.na(result_zeros$R)))

  # Single case
  single_case <- c(1, rep(0, length(test_data$dates) - 1))
  result_single <- wallinga_lipsitch(single_case, test_data$dates, test_data$si_mean, test_data$si_sd)
  # Should have R = NA for the single case (no secondary cases to estimate R)
  expect_true(is.na(result_single$R[1]))
})
