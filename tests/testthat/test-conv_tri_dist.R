test_that("check route argument", {
  expect_error(conv_tri_dist(x = 1:10, sigma = 3, r = 10, mu = 5, route = 0))
  expect_error(conv_tri_dist(x = 1:10, sigma = 3, r = 10, mu = 5, route = ""))
})

test_that("check quantity argument", {
  expect_error(conv_tri_dist(x = 1:10, sigma = 3, r = 10, mu = 5, route = 1, quantity = 1))
  expect_error(conv_tri_dist(x = 1:10, sigma = 3, r = 10, mu = 5, route = 1, quantity = ""))
})
