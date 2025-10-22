get_sill <- function(model) {
  if(methods::is(model, "variogramModel"))
    gstat::variogramLine(model, dist_vector = 0, covariance = TRUE)$gamma
  else if(methods::is(model, "StVariogramModel"))
    gstat::variogramSurface(model, data.frame(spacelag = 0, timelag = 0), covariance = TRUE)$gamma
}

#' Spatial Ordinary or Universal Kriging
#'
#' @param formula formula that defines the response variable as a linear model of covariates
#' @param data object of class `Spatial` or `sf` that contains the response variable, covariates, and coordinates
#' @param newdata object of class `Spatial` or `sf` with prediction locations; should contain attributes with the covariates (if any)
#' @param model variogram model
#' @param errorVar variances of measurement errors
#' @param computeVar logical; compute prediction variances or not
#' @param epsilon small value added to nugget to avoid singular covariance matrix
#'
#' @references Hsu CC, Tsai DR, Su SY, Jhuang JR, Chiang CJ, Yang YW, Lee WC.
#'     A Stabilized Kriging Method for Mapping Disease Rates.
#'     Journal of Epidemiol (2023); 33(4):201-208.
#'
#' @export
#'
#' @importFrom methods as
#'
krige_sp <- function(formula, data, newdata, model, errorVar, computeVar = FALSE, epsilon = 1e-7) {

  mf <- match.call()
  ind <- match(c("formula", "data", "errorVar"), names(mf), 0L)
  mf <- mf[c(1L, ind)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  if(inherits(data, c("sf", "sfc"))) data <- as(data, "Spatial")
  if(inherits(newdata, c("sf", "sfc"))) newdata <- as(newdata, "Spatial")

  dmat <- stats::model.matrix(formula, mf)
  dmat.new <- stats::model.matrix(stats::update(formula, NULL ~ .), newdata)

  nug <- epsilon
  nug_ind <- which(model$model %in% c("Nug", "Err"))
  if(length(nug_ind)) {
    nug <- pmax(model$psill[nug_ind], epsilon)
    model <- model[-nug_ind, ]
  }

  if(model$psill == 0) model <- gstat::vgm(epsilon, "Nug", 0)

  covar <- gstat::variogramLine(model, dist_vector = sp::spDists(data), covariance = TRUE)
  n <- nrow(data@data)

  var <- stats::model.extract(mf, "errorVar")
  if(is.null(var)) var <- nug

  C <- rbind(
    cbind(covar + diag(var, n), dmat),
    cbind(t(dmat), diag(0, ncol(dmat)))
  )

  C.inv <- solve(C)

  covar0 <- gstat::variogramLine(model, dist_vector = sp::spDists(data, newdata), covariance = TRUE)
  D <- rbind(covar0, t(dmat.new))
  w <- C.inv %*% D
  if(!all(abs(colSums(w[1:n, , drop = FALSE]) - 1) < 1e-5))
    warning("The kriging weights do not sum to 1")

  class.new <- sub("DataFrame$", "", class(newdata))
  newdata.2 <- as(newdata, class.new)
  pred <- data.frame(var1.pred = c(t(stats::model.response(mf)) %*% w[1:n, , drop = FALSE]), row.names = row.names(newdata))
  result <- sp::SpatialPointsDataFrame(newdata.2, pred)
  result <- as(result, paste0(class.new, "DataFrame"))

  if(grepl("Point|Pixel", class.new))
    result@coords.nrs <- seq_len(nrow(newdata@bbox))

  if(computeVar) {
    result$var1.var <- get_sill(model) - colSums(w * D)
  }

  return(result)
}
