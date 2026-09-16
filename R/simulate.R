dist_ellipse <- function(init_dist, direction, anis) {
  d1 <- init_dist * cos((direction - anis[1]) * (pi/180))
  d2 <- init_dist * sin((direction - anis[1]) * (pi/180))
  sqrt(d1^2 + (d2 / anis[2])^2)
}

# map angles (radian) to [-pi, pi)
wrap_angle <- function(angle) {
  atan2(sin(angle), cos(angle))
}

shift_center <- function(c1, c2, t, t_lim) {
  c1.coord <- sp::coordinates(c1)
  c2.coord <- sp::coordinates(c2)
  c.coord <- t(c(c1.coord) + outer(c(c2.coord - c1.coord), t - min(t_lim)) / (max(t_lim) - min(t_lim)))
  c <- sp::SpatialPoints(c.coord, proj4string = c1@proj4string)
  return(c)
}

gauss2d <- function(s, c0, bw, anis) {
  theta <- geosphere::bearing(s, c0)
  d <- sp::spDists(s, c0)
  dd <- dist_ellipse(d, theta, anis)
  exp(-dd^2 / (bw / anis[2]))
}

#' Spatiotemporal Random Surface Generator
#'
#' @param sgrid spatial grid
#' @param tgrid temporal grid
#' @param vgrid grid of values multiplied to the kernels
#' @param anis_angle anisotropy angle, the angle of the main axis direction
#' @param anis_ratio anisotropy ratio, the ratio of the minor range to the major range
#' @param sbw spatial bandwidth
#' @param tbw temporal bandwidth
#' @param n number of grid points to sample
#' @param seed random seed
#'
#' @return a function that takes spatial and temporal coordinates and returns values of the random surface
#'
#' @export
#'
st_rand_surface <- function(sgrid, tgrid, vgrid = c(-1, 1),
                            anis_angle = 1:180, anis_ratio = 1,
                            sbw, tbw, n = 101, seed) {

  if(!is(sgrid, "Spatial")) sgrid <- as(sgrid, "Spatial")
  if(!missing(seed)) set.seed(seed)

  sx <- sample(length(sgrid), n, replace = TRUE)
  tx <- sample(        tgrid, n, replace = TRUE)
  vx <- sample(        vgrid, n, replace = TRUE)
  ax <- sample(   anis_angle, n, replace = TRUE)
  rx <- sample(   anis_ratio, n, replace = TRUE)

  function(s, t) {
    if(!is(s, "Spatial")) s <- as(s, "Spatial")

    Map(\(si, ti, vi, ai, ri) {
      ss <- gauss2d(s, sgrid[si, ], bw = sbw, anis = c(ai, ri))
      tt <- exp(-(t - ti)^2 / tbw)
      vi * outer(c(ss), tt, FUN = "*")
    }, sx, tx, vx, ax, rx) |> Reduce(f = `+`)
  }
}

#' Spatiotemporal Risk Simulation
#'
#' @param s spatial coordinate
#' @param t index of time groups (a vector)
#' @param age index of age groups (a single value)
#' @param c1 center of the first static hotspot
#' @param c2 center of the second static hotspot
#' @param c_start starting center of the dynamic hotspot
#' @param c_end ending center of the dynamic hotspot
#' @param bw_stat bandwidth of the static hotspot(s)
#' @param bw_dyn bandwidth of the dynamic hotspot
#' @param anis_stat anisotropy parameters of the static hotspot(s); see notes of [gstat::vgm()]
#' @param anis_dyn anisotropy parameters of the dynamic hotspot; see notes of [gstat::vgm()]
#' @param rr_stat ratio of the highest and lowest disease rates for the static hotspot(s)
#' @param rr_dyn ratio of the highest and lowest disease rates for the dynamic hotspot
#' @param t_lim time limits
#' @param t_turn turning point in time
#' @param t_par time parameter
#' @param n_arm number of spiral arms
#' @param arm_bw bandwidth of spiral arms
#' @param arm_tight spiral tightness parameter
#' @param arm_eye size of the spiral center as a fraction of `bw_dyn`
#' @param arm_rot angular velocity of the spinning spiral (in degrees per unit time)
#' @param age_fun a function of age that outputs age-specific baseline rates (in log scale)
#' @param scenario simulation scenario
#'   1. convergence–dissipation
#'   2. outward rippling
#'   3. expansion–contraction gaussian
#'   4. expansion–contraction spinning spiral
#' @param ... arguments to be passed to \code{st_risk()}
#'
#' @return a matrix of disease rates in log scale
#'
#' @name simulation
#' @export
#'
st_risk <- function(s, t, age,
                    c1, c2 = c1, c_start, c_end,
                    bw_stat, bw_dyn,
                    anis_stat = c(0, 1), anis_dyn = c(0, 1),
                    rr_stat, rr_dyn,
                    t_lim, t_turn, t_par,
                    n_arm, arm_bw, arm_tight, arm_eye, arm_rot,
                    age_fun, scenario) {

  if(!is(s, "Spatial")) s <- as(s, "Spatial")

  bw_stat <- rep_len(bw_stat, 2)
  rr_stat <- rep_len(rr_stat, 2)
  if(!is.list(anis_stat)) anis_stat <- list(anis_stat)
  anis_stat <- rep_len(anis_stat, 2)

  int1 <- log(rr_stat[1]) * gauss2d(s, c1, bw = bw_stat[1], anis = anis_stat[[1]])
  int2 <- log(rr_stat[2]) * gauss2d(s, c2, bw = bw_stat[2], anis = anis_stat[[2]])
  int <- (int1 + int2) / 2

  center <- shift_center(c_start, c_end, t, t_lim)
  theta <- sapply(seq_along(center), \(i) geosphere::bearing(s, center[i]))
  d <- sp::spDists(s, center)
  static <- age_fun(age) + array(int, dim(d))
  dd <- dist_ellipse(d, theta, anis_dyn)

  if(scenario %in% 1:3) {
    xx <- switch(scenario, {
      # 1
      exp(-(t(dd)^2 / (bw_dyn / anis_dyn[2]) + (t - t_turn)^2 / t_par))
    }, {
      # 2
      r <- t_par * pmax(0, t - t_turn)
      exp(-((t(dd) - r)^2 / (bw_dyn / anis_dyn[2])))
    }, {
      # 3
      h <- bw_dyn * exp(-(t - t_turn)^2 / t_par)
      exp(-t(dd)^2 / (h / anis_dyn[2]))
    })

    res <- static + t(log(rr_dyn) * xx)
  }

  if(scenario == 4L) {
    n_arm <- as.integer(n_arm)
    res <- 0

    for (k in 0:(n_arm - 1)) {
      dtheta <- t(theta*(pi/180) - arm_tight*log(dd) - 2*pi*k/n_arm)
      res <- res + exp(-(wrap_angle(dtheta - arm_rot*(pi/180) * t) / arm_bw)^2)
    }

    h <- bw_dyn * exp(-(t - t_turn)^2 / t_par)
    eye <- exp(-dd^2 / (bw_dyn * arm_eye))
    res <- pmax(res, t(eye))
    res <- static + t(res * log(rr_dyn) * exp(-(t(dd)^2 / (h / anis_dyn[2]))))
  }

  return(res)
}

#' @rdname simulation
#' @export
#'
st_risk1 <- function(...) {
  st_risk(..., scenario = 1L)
}

#' @rdname simulation
#' @export
#'
st_risk2 <- function(...) {
  st_risk(..., scenario = 2L)
}

#' @rdname simulation
#' @export
#'
st_risk3 <- function(...) {
  st_risk(..., scenario = 3L)
}

#' @rdname simulation
#' @export
#'
st_risk4 <- function(..., n_arm = 3L, arm_bw = 0.5, arm_tight = 1.0, arm_eye = 0.01, arm_rot = 0) {
  st_risk(...,
          n_arm = n_arm,
          arm_bw = arm_bw,
          arm_tight = arm_tight,
          arm_eye = arm_eye,
          arm_rot = arm_rot,
          scenario = 4L)
}
