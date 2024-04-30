test_that("check function arguments", {
  expect_error(weighted_var("a", 1:10))
  expect_error(weighted_var(1, "b"))
})
