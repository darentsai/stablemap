#' Symmetric Mean Absolute Percentage Error
#'
#' @param x predicted values
#' @param y true values
#'
#' @return A numeric value
#' @export
#'
smape <- function(x, y) {
  mean(abs(x - y) / (x + y), na.rm = TRUE) * 100
}

#' Root Mean Square Error
#'
#' @param x predicted values
#' @param y true values
#'
#' @return A numeric value
#' @export
#'
rmse <- function(x, y) {
  sqrt(mean((x - y)^2, na.rm = TRUE))
}
