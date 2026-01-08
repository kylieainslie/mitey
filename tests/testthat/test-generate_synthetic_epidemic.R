# Test file for generate_synthetic_epidemic function

test_that("generate_synthetic_epidemic returns correct structure", {
  set.seed(123)
  true_r <- rep(1.5, 30)

  result <- generate_synthetic_epidemic(
    true_r = true_r,
    si_mean = 5,
    si_sd = 2,
    si_dist = "gamma"
  )

  # Check return type and structure

  expect_s3_class(result, "data.frame")
  expect_named(result, c("date", "true_r", "incidence"))
  expect_equal(nrow(result), 30)

  # Check column types
  expect_s3_class(result$date, "Date")
  expect_type(result$true_r, "double")
  expect_true(is.numeric(result$incidence))

  # Check that true_r is preserved
  expect_equal(result$true_r, true_r)

  # Check that dates are sequential
  expect_equal(as.numeric(diff(result$date)), rep(1, 29))
})

test_that("generate_synthetic_epidemic works with different distributions", {
  set.seed(456)
  true_r <- rep(1.2, 20)

  # Test gamma distribution
  result_gamma <- generate_synthetic_epidemic(
    true_r = true_r,
    si_mean = 7,
    si_sd = 3,
    si_dist = "gamma"
  )
  expect_equal(nrow(result_gamma), 20)
  expect_true(all(result_gamma$incidence >= 0))

  # Test normal distribution
  result_normal <- generate_synthetic_epidemic(
    true_r = true_r,
    si_mean = 7,
    si_sd = 3,
    si_dist = "normal"
  )
  expect_equal(nrow(result_normal), 20)
  expect_true(all(result_normal$incidence >= 0))
})

test_that("generate_synthetic_epidemic respects initial_cases parameter", {
  set.seed(789)
  true_r <- rep(1.0, 10)

  result <- generate_synthetic_epidemic(
    true_r = true_r,
    si_mean = 5,
    si_sd = 2,
    si_dist = "gamma",
    initial_cases = 50
  )

  # First day should have initial_cases
  expect_equal(result$incidence[1], 50)
})

test_that("generate_synthetic_epidemic produces reasonable epidemics", {
  set.seed(111)

  # Growing epidemic (R > 1)
  growing_r <- rep(2.0, 50)
  growing_epidemic <- generate_synthetic_epidemic(
    true_r = growing_r,
    si_mean = 5,
    si_sd = 2,
    si_dist = "gamma",
    initial_cases = 10
  )

  # Total cases should generally increase over time
  early_cases <- sum(growing_epidemic$incidence[1:10])
  late_cases <- sum(growing_epidemic$incidence[41:50])
  expect_true(late_cases > early_cases)

  # Declining epidemic (R < 1)
  set.seed(222)
  declining_r <- rep(0.5, 50)
  declining_epidemic <- generate_synthetic_epidemic(
    true_r = declining_r,
    si_mean = 5,
    si_sd = 2,
    si_dist = "gamma",
    initial_cases = 100
  )

  # Cases should generally decrease
  early_declining <- sum(declining_epidemic$incidence[1:10])
  late_declining <- sum(declining_epidemic$incidence[41:50])
  expect_true(late_declining <= early_declining)
})

test_that("generate_synthetic_epidemic handles edge cases", {
  set.seed(333)

  # Few days
  few_days <- generate_synthetic_epidemic(
    true_r = rep(1.5, 5),
    si_mean = 5,
    si_sd = 2,
    si_dist = "gamma"
  )
  expect_equal(nrow(few_days), 5)
  expect_equal(few_days$incidence[1], 10)  # Default initial_cases

  # Very short serial interval
  short_si <- generate_synthetic_epidemic(
    true_r = rep(1.2, 20),
    si_mean = 2,
    si_sd = 0.5,
    si_dist = "gamma",
    initial_cases = 10
  )
  expect_equal(nrow(short_si), 20)
  expect_true(all(short_si$incidence >= 0))

  # Very long serial interval
  long_si <- generate_synthetic_epidemic(
    true_r = rep(1.2, 20),
    si_mean = 14,
    si_sd = 5,
    si_dist = "gamma",
    initial_cases = 10
  )
  expect_equal(nrow(long_si), 20)
  expect_true(all(long_si$incidence >= 0))
})

test_that("generate_synthetic_epidemic is reproducible with seed", {
  set.seed(999)
  result1 <- generate_synthetic_epidemic(
    true_r = rep(1.5, 30),
    si_mean = 5,
    si_sd = 2,
    si_dist = "gamma"
  )

  set.seed(999)
  result2 <- generate_synthetic_epidemic(
    true_r = rep(1.5, 30),
    si_mean = 5,
    si_sd = 2,
    si_dist = "gamma"
  )

  expect_equal(result1$incidence, result2$incidence)
})
