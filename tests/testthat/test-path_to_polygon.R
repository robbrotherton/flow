test_that("path_to_polygon works", {

  path <- data.frame(x = c(1, 2, 3, 4, 5),
                     y = c(0, 0, 0, 0, 0),
                     a = c(0, 0, 0, 0, 0))

  expected_output <- data.frame(x = c(1:5, 5:1),
                                y = c(rep(1, 5), rep(-1, 5)))

  expect_equal(path_to_polygon(path), expected_output)

})
