# Test cases for integrate_all_components function
test_that("integrate_components_wrapper calculates integrals correctly", {
  # Test case 1: Data point = 0
  result1 <- integrate_components_wrapper(0, 10, 2)
  expect_length(result1, 7)
  expect_true(all(result1 >= 0))  # Check that all values are non-negative

  # Test case 2: Data point = 1
  result2 <- integrate_components_wrapper(1, 10, 2)
  expect_length(result2, 7)
  expect_true(all(result2 >= 0))  # Check that all values are non-negative

  # Add more test cases as needed
})
