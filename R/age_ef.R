#' Age Effect
#'
#' @param n_age number of age groups
#' @param min_rate rate in the first age group
#' @param rr ratio of the rate in the last age group to that in the first age group
#'
#' @return a function that takes integer indices of age groups and returns the logarithmic rates
#'
#' @name age_ef
#' @export
#'
exp_age_ef <- function(n_age, min_rate = 25e-5, rr = 10) {
  function(age) {
    if(n_age == 1)
      log(min_rate)
    else
      log(min_rate) + log(rr) / (n_age - 1) * (age - 1)
  }
}
