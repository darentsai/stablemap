path <- function(c1, c2, t, t_n) {
  c1.coord <- sp::coordinates(c1)
  c2.coord <- sp::coordinates(c2)
  c.coord <- t(c(c1.coord) + outer(c(c2.coord - c1.coord), t - 1) / (t_n - 1))
  c <- sp::SpatialPoints(c.coord, proj4string = c1@proj4string)
  return(c)
}

static_risk <- function(s, c1, c2, bw) {
  d <- sp::spDists(s, c1)
  logRR <- exp(-d^2 / bw)
  if(!identical(c1, c2)) {
    d <- sp::spDists(s, c2)
    logRR <- (logRR + exp(-d^2 / bw)) / 2
  }
  return(log(2) * logRR)
}

#' Exponential Age Effect
#'
#' @param age index of age groups
#' @param age_n number of age groups
#' @param min_rate rate in the first age group
#' @param ratio ratio of the rate in the last age group to that in the first age group
#'
#' @export
#'
exp_age_ef <- function(age, age_n = 9, min_rate = 25e-5, ratio = 10) {
  log(min_rate) + log(ratio) / (age_n - 1) * (age - 1)
}

#' Spatiotemporal Simulation
#'
#' @param s spatial coordinate
#' @param t index of time groups
#' @param age index of age groups
#' @param c1 center of the first static hotspot
#' @param c2 center of the second static hotspot
#' @param c_start starting center of the dynamic hotspot
#' @param c_end ending center of the dynamic hotspot
#' @param s_bw_static bandwidth of static hotspot
#' @param s_bw_dynamic bandwidth of dynamic hotspot
#' @param t_n number of time groups
#' @param t_turn turning time point
#' @param t_bw bandwidth in time
#' @param age_fun age effect function
#'
#' @export
#' @name simulation
#'
scenario1 <- function(s, t, age,
                      c1, c2, c_start, c_end,
                      s_bw_static, s_bw_dynamic,
                      t_n, t_turn, t_bw,
                      age_fun = exp_age_ef) {
  int <- static_risk(s, c1, c2, s_bw_static)
  d <- sp::spDists(s, path(c_start, c_end, t, t_n))
  age_fun(age) + array(int, dim(d)) + t(exp(-(t - t_turn)^2 / t_bw) * exp(-t(d)^2 / s_bw_dynamic))
}

#' @export
#' @rdname simulation
#'
scenario2 <- function(s, t, age,
                      c1, c2, c_start, c_end,
                      s_bw_static, s_bw_dynamic,
                      t_n, t_turn, t_bw,
                      age_fun = exp_age_ef) {
  int <- static_risk(s, c1, c2, s_bw_static)
  d <- sp::spDists(s, path(c_start, c_end, t, t_n))
  r <- t_bw * pmax(0, t - t_turn)
  age_fun(age) + array(int, dim(d)) + t(exp(-(t(d) - r)^2 / s_bw_dynamic))
}

#' @export
#' @rdname simulation
#'
scenario3 <- function(s, t, age,
                      c1, c2, c_start, c_end,
                      s_bw_static, s_bw_dynamic,
                      t_n, t_turn, t_bw,
                      age_fun = exp_age_ef) {
  int <- static_risk(s, c1, c2, s_bw_static)
  d <- sp::spDists(s, path(c_start, c_end, t, t_n))
  h <- s_bw_dynamic * exp(-(t - t_turn)^2 / t_bw)
  age_fun(age) + array(int, dim(d)) + t(exp(-t(d)^2 / h))
}
