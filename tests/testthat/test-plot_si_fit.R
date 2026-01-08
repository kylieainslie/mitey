# Test file for plot_si_fit function

test_that("plot_si_fit returns a ggplot object for normal distribution", {
  set.seed(123)
  icc_data <- c(
    round(rnorm(20, mean = 2, sd = 1)),
    round(rnorm(50, mean = 12, sd = 3)),
    round(rnorm(20, mean = 24, sd = 4))
  )
  icc_data <- pmax(icc_data, 0)

  p <- plot_si_fit(
    dat = icc_data,
    mean = 12.5,
    sd = 3.2,
    weights = c(0.2, 0.6, 0.15, 0.05),
    dist = "normal"
  )

  expect_s3_class(p, "ggplot")

  # Check that plot has expected layers
  expect_true(length(p$layers) >= 2)  # At least histogram and stat_function
})

test_that("plot_si_fit returns a ggplot object for gamma distribution", {
  set.seed(456)
  icc_data <- round(rgamma(100, shape = 16, rate = 2))

  p <- plot_si_fit(
    dat = icc_data,
    mean = 8.0,
    sd = 2.0,
    weights = c(0.25, 0.65, 0.10),
    dist = "gamma"
  )

  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 2)
})

test_that("plot_si_fit handles scaling_factor parameter", {
  set.seed(789)
  icc_data <- round(rnorm(50, mean = 10, sd = 3))
  icc_data <- pmax(icc_data, 1)

  # Default scaling
  p1 <- plot_si_fit(
    dat = icc_data,
    mean = 10,
    sd = 3,
    weights = c(0.1, 0.7, 0.15, 0.05),
    dist = "normal",
    scaling_factor = 1
  )

  # Custom scaling
  p2 <- plot_si_fit(
    dat = icc_data,
    mean = 10,
    sd = 3,
    weights = c(0.1, 0.7, 0.15, 0.05),
    dist = "normal",
    scaling_factor = 0.5
  )

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_si_fit works with real si_estim output", {
  # Simulate data similar to what would be used with si_estim
  set.seed(111)
  icc_intervals <- c(
    rep(6, 4), rep(7, 8), rep(8, 14), rep(9, 31),
    rep(10, 29), rep(11, 42), rep(12, 25), rep(13, 16),
    rep(14, 16), rep(15, 10), rep(16, 4), rep(17, 2), rep(18, 2)
  )

  # Run si_estim to get real parameters
  si_results <- si_estim(icc_intervals)

  # Create plot with real results
  p <- plot_si_fit(
    dat = icc_intervals,
    mean = si_results$mean,
    sd = si_results$sd,
    weights = c(
      si_results$wts[1],
      si_results$wts[2] + si_results$wts[3],
      si_results$wts[4] + si_results$wts[5],
      si_results$wts[6] + si_results$wts[7]
    ),
    dist = "normal"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_si_fit handles small datasets", {
  small_data <- c(5, 6, 7, 8, 9)

  p <- plot_si_fit(
    dat = small_data,
    mean = 7,
    sd = 1.5,
    weights = c(0.1, 0.8, 0.08, 0.02),
    dist = "normal"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_si_fit can be built without errors", {
  set.seed(222)
  icc_data <- round(rnorm(50, mean = 10, sd = 3))
  icc_data <- pmax(icc_data, 1)

  p <- plot_si_fit(
    dat = icc_data,
    mean = 10,
    sd = 3,
    weights = c(0.1, 0.7, 0.15, 0.05),
    dist = "normal"
  )

  # Check that the plot can be built (this would catch issues with stat_function)
  built <- ggplot2::ggplot_build(p)
  expect_type(built, "list")
  expect_true("data" %in% names(built))
})
