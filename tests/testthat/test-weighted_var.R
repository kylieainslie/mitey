test_that("check function arguments", {
  expect_error(weighted_var("a", 1:10))
  expect_error(weighted_var(1, "b"))
  expect_error(weighted_var(1:11, seq(0.1, 1, length.out = 10)))
})

test_that("check output", {
  expect_equal(weighted_var(1:10, seq(0.1, 1, length.out = 10)), 6.875)
  expect_equal(weighted_var(1,1), NaN)
})
