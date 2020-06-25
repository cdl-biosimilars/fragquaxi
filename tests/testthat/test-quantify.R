test_that("rt limits are correctly identified", {
  ms_header <- tibble::tibble(
    seqNum = 1:7,
    retentionTime = c(2, 4, 7, 7.5, 8, 11, 30),
    msLevel = c(1L, 2L, 1L, 1L, 2L, 1L, 1L)
  )
  rt_limits <- tibble::tibble(
    rt_min = c(-Inf, -Inf, 3.9, 9, 3),
    rt_max = c(+Inf, 8, +Inf, 10, 15)
  )
  expect_equal(
    find_rt_limits(ms_header, rt_limits),
    tibble::tibble(
      rt_min = c(2, 2, 3.9, 9, 3),
      rt_max = c(30, 8, 30, 10, 15),
      scans = list(1:7, 1:5, 2:7, NA_integer_, 2:6)
    )
  )
  expect_equal(
    find_rt_limits(ms_header, c(3, 7)),
    tibble::tibble(rt_min = 3, rt_max = 7, scans = list(2:3))
  )
  expect_equal(
    find_rt_limits(ms_header[ms_header$msLevel == 1,], rt_limits),
    tibble::tibble(
      rt_min = c(2, 2, 3.9, 9, 3),
      rt_max = c(30, 8, 30, 10, 15),
      scans = list(c(1, 3, 4, 6, 7), c(1, 3, 4), c(3, 4, 6, 7),
                   NA_integer_, c(3, 4, 6))
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
