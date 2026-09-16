#' Symmetric Mean Absolute Percentage Error (SMAPE)
#'
#' Calculate the symmetric mean absolute percentage error (SMAPE).
#'
#' @param y_true true values
#' @param y_pred predicted values
#' @param na.rm logical; whether NA values should be removed before the computation proceeds
#'
#' @return A numeric value
#'
#' @export
#'
smape <- function(y_true, y_pred, na.rm = FALSE) {
  mean(abs(y_pred - y_true) / (abs(y_pred) + abs(y_true)), na.rm = na.rm)
}

#' Poisson Deviance-based \eqn{R^2} for Count Data
#'
#' Calculate the Poisson deviance-based pseudo-\eqn{R^2} for count data.
#'
#' @param y_pred predicted values
#' @param y_true true values
#' @param na.rm logical; whether NA values should be removed before the computation proceeds
#'
#' @references A. Colin Cameron and Frank A. G. Windmeijer.
#'     R-Squared Measures for Count Data Regression Models With Applications to Health-Care Utilization.
#'     Journal of Business & Economic Statistics (1996); 14(2):209-220.
#'
#' @return A numeric value
#'
#' @export
#'
r2pois <- function(y_true, y_pred, na.rm = FALSE) {
  y_bar <- mean(y_true, na.rm = na.rm)
  dev0 <- devPois(y_true, y_bar, na.rm = na.rm)
  dev1 <- devPois(y_true, y_pred, na.rm = na.rm)
  1 - dev1 / dev0
}

devPois <- function(y_true, y_pred, na.rm = FALSE) {
  sum(ifelse(y_true == 0, 0, y_true * log(y_true / y_pred)) - (y_true - y_pred), na.rm = na.rm)
}
