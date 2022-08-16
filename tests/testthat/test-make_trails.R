test_that("make_trails works", {
  ff <- make_flowfield(angle = 0)

  trail <- make_trails(ff,
                       data.frame(x = 50, y = 50))

  expect_equal(trail,
               data.frame(x = 1:99,
                          y = rep(50, 99),
                          group = 0))
})
