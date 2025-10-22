path <- function(c1, c2, t, t_n) {
  c1.coord <- sp::coordinates(c1)
  c2.coord <- sp::coordinates(c2)
  c.coord <- t(c(c1.coord) + outer(c(c2.coord - c1.coord), t - 1) / (t_n - 1))
  c <- sp::SpatialPoints(c.coord, proj4string = c1@proj4string)
  return(c)
}

static_risk <- function(s, c1, c2, bw, rr) {
  d <- sp::spDists(s, c1)
  res <- exp(-d^2 / bw)
  if(!identical(c1, c2)) {
    d <- sp::spDists(s, c2)
    res <- (res + exp(-d^2 / bw)) / 2
  }
  return(log(rr) * res)
}

#' Exponential Age Effect
#'
#' @param age_n number of age groups
#' @param min_rate rate in the first age group
#' @param rr ratio of the rate in the last age group to that in the first age group
#'
#' @export
#' @name age_ef
#'
exp_age_ef <- function(age_n, min_rate = 25e-5, rr = 10) {
  function(age) {
    log(min_rate) + log(rr) / (age_n - 1) * (age - 1)
  }
}

#' Spatiotemporal Simulation
#'
#' @param s spatial coordinate
#' @param t index of time groups (a vector)
#' @param age index of age groups (a single value)
#' @param c1 center of the first static hotspot
#' @param c2 center of the second static hotspot
#' @param c_start starting center of a dynamic hotspot
#' @param c_end ending center of a dynamic hotspot
#' @param s_bw_static bandwidth of a static hotspot
#' @param s_bw_dynamic bandwidth of a dynamic hotspot
#' @param rr_static ratio of the highest and lowest disease rates for a static hotspot
#' @param rr_dynamic ratio of the highest and lowest disease rates for a dynamic hotspot
#' @param t_n number of time groups
#' @param t_turn turning time point
#' @param t_bw bandwidth in time
#' @param age_fun a function of age that outputs age-specific baseline rates
#' @param scenario simulation scenario
#' \itemize{
#'   \item 1: convergence–dissipation
#'   \item 2: outward rippling
#'   \item 3: expansion–contraction
#' }
#' @param ... arguments to be passed to \code{st_sim()}
#'
#' @return A matrix of logarithmic disease rates
#'
#' @export
#' @name simulation
#'
st_sim <- function(s, t, age,
                   c1, c2, c_start, c_end,
                   s_bw_static, s_bw_dynamic,
                   rr_static, rr_dynamic,
                   t_n, t_turn, t_bw,
                   age_fun, scenario) {

  if(inherits(s, c("sf", "sfc"))) s <- as(s, "Spatial")

  int <- static_risk(s, c1, c2, s_bw_static, rr_static)
  d <- sp::spDists(s, path(c_start, c_end, t, t_n))
  static <- age_fun(age) + array(int, dim(d))

  if(scenario == 1L) {
    res <- static + t(log(rr_dynamic) * exp(-(t - t_turn)^2 / t_bw - t(d)^2 / s_bw_dynamic))
  } else if(scenario == 2L) {
    r <- t_bw * pmax(0, t - t_turn)
    res <- static + t(log(rr_dynamic) * exp(-(t(d) - r)^2 / s_bw_dynamic))
  } else if(scenario == 3L) {
    h <- s_bw_dynamic * exp(-(t - t_turn)^2 / t_bw)
    res <- static + t(log(rr_dynamic) * exp(-t(d)^2 / h))
  }

  return(res)
}

#' @export
#' @rdname simulation
#'
st_sim1 <- function(...) {
  st_sim(..., scenario = 1L)
}

#' @export
#' @rdname simulation
#'
st_sim2 <- function(...) {
  st_sim(..., scenario = 2L)
}

#' @export
#' @rdname simulation
#'
st_sim3 <- function(...) {
  st_sim(..., scenario = 3L)
}
