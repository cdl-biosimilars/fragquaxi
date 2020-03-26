#' Trapezoidal integration with limits.
#'
#' Compute the area of a function with values y at the points x within given
#' limits.
#'
#' @param x x-coordinates of points on the x-axis.
#' @param y Function values associated with the x-coordinates.
#' @param lower Lower integration boundary.
#' @param upper Upper integration boundary.
#'
#' @return Area under the curve defined by the (x, y)-tuples.
#' @export
#'
#' @examples
#' x <- 0:5
#' y <- c(1, 4, 2, 3, 3, 0)
#'
#' trapz_limits(x, y, -2, 10)
#'
#' trapz_limits(x, y, -2, -1)
#'
#' trapz_limits(x, y, -2, .4)
#'
#' trapz_limits(x, y, .4, .7)
#'
#' trapz_limits(x, y, .4, 3.6)
#'
#' trapz_limits(x, y, .4, 7)
#'
#' trapz_limits(x, y, 6, 7)
#'
#' \dontrun{
#' # throws an error
#' trapz_limits(x, y, 4, 3)
#' }
trapz_limits <- function(x, y, lower, upper) {
  if (upper < lower)
    stop("'upper' must be larger than 'lower'")

  y_lim <- stats::approx(x, y, c(lower, upper))$y
  ind <- which(dplyr::between(x, lower, upper))
  x_min <- x[1]
  x_max <- x[length(x)]

  if (lower > x_max | upper < x_min)
    return(0)

  if (lower > x_min) {
    if (upper < x_max) {
      x <- c(lower, x[ind], upper)
      y <- c(y_lim[1], y[ind], y_lim[2])
    } else {
      x <- c(lower, x[ind])
      y <- c(y_lim[1], y[ind])
    }
  } else {
    if (upper < x_max) {
      x <- c(x[ind], upper)
      y <- c(y[ind], y_lim[2])
    }
  }

  sum(
    (dplyr::lead(x) - x) * (dplyr::lead(y) + y) / 2,
    na.rm = TRUE
  )
}
