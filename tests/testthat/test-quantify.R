test_that("rt limits are correctly identified", {
  rt <- c(2, 4, 7, 7.5, 8, 11, 30)
  rt_limits <- tibble::tribble(
    ~rt_min, ~rt_max,
    -Inf, +Inf,
    -Inf, 8,
    3.9, +Inf,
    9, 10,
    3, 15,
  )
  expect_equal(
    find_rt_limits(rt, rt_limits),
    tibble::tribble(
      ~rt_min, ~rt_max, ~scan_min,   ~scan_max,
      2,       30,      1L,          7L,
      2,        8,      1L,          5L,
      3.9,     30,      2L,          7L,
      9,       10,      NA_integer_, NA_integer_,
      3,       15,      2L,          6L
    )
  )
  expect_equal(
    find_rt_limits(rt, c(3, 7)),
    tibble::tribble(
      ~rt_min, ~rt_max, ~scan_min,   ~scan_max,
      3,       7,       2L,          3L,
    )
  )
})

test_that("trapezoidal integration works", {
  trapz_limits_p <- purrr::partial(
    trapz_limits, x = 0:5, y = c(1, 4, 2, 3, 3, 0)
  )
  expect_equal(trapz_limits_p(-2, 10), 12.5)
  expect_equal(trapz_limits_p(-2, -1), 0)
  expect_equal(trapz_limits_p(-2, .4), 0.64)
  expect_equal(trapz_limits_p(.4, .7), 0.795)
  expect_equal(trapz_limits_p(.4, 3.6), 9.16)
  expect_equal(trapz_limits_p(.4, 7), 11.86)
  expect_equal(trapz_limits_p(6, 7), 0)
  expect_error(trapz_limits_p(4, 3))
})
