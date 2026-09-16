near0 <- 1e-32

get_vgm_param <- function(model, covariance) {
  if(is(model, "variogramModel")) {
    gstat::variogramLine(model, covariance = covariance, dist_vector = near0)$gamma
  } else if(is(model, "StVariogramModel")) {
    gstat::variogramSurface(
      model, covariance = covariance,
      dist_grid = data.frame(spacelag = near0, timelag = near0)
    )$gamma
  }
}

#' Extract Nugget and Partial Sill from a Variogram Model
#'
#' @param model object of class `variogramModel` or `StVariogramModel`
#'
#' @export
#'
get_nugget <- function(model) {
  get_vgm_param(model, FALSE)
}

#' @rdname get_nugget
#' @export
#'
get_sill <- function(model) {
  get_vgm_param(model, TRUE)
}

#' Slice Plot of Space-Time Variogram Models
#'
#' @param x object of class `StVariogram`
#' @param model object of class `StVariogramModel`
#' @param aspect slicing aspect
#' @param col line colors
#' @param ... arguments passed to [ggplot2::facet_wrap()]
#'
#' @export
#'
slice_stvgm <- function(x, model = NULL, aspect = 1, col = sp::bpy.colors(), ...) {

  if(is(model, "StVariogramModel"))
    model <- list(model)
  if(!is.null(model)) {
    mod.name <- sapply(model, \(m) m$stModel)
    names(model) <- tail(make.names(c(unique(mod.name[duplicated(mod.name)]), mod.name), unique = TRUE),
                         length(mod.name))
  }

  asp <- c("spacelag", "timelag")[seq(aspect, -aspect+3)]
  x$spacelag <- x$avgDist
  x <- x[-1, c(asp, "gamma")]

  x$flag <- FALSE
  jump <- data.frame(1e-5, unique(x[[asp[2]]]), flag = TRUE)
  names(jump)[1:2] <- asp
  dat <- dplyr::bind_rows(c(
    list(sample = x),
    lapply(model, gstat::variogramSurface, dist_grid = rbind(jump, x[-3]))
  ), .id = "model")
  dat$model <- factor(dat$model, levels = unique(dat$model))

  ggplot2::ggplot(dat, aes(x = .data[[asp[1]]], y = gamma, colour = .data[[asp[2]]], group = .data[[asp[2]]])) +
    ggplot2::geom_line() +
    ggplot2::geom_point(shape = 1, data = dat[!dat$flag, ]) +
    ggplot2::scale_color_gradientn(colours = col) +
    ggplot2::facet_wrap(~ model, ...) +
    ggplot2::labs(y = expression(paste("semi-variogram ", (gamma)))) +
    ggplot2::theme_bw()
}
