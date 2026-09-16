#' Spatial Kriging (Ordinary or Universal)
#'
#' @param formula formula that defines the response variable as a linear model of covariates
#' @param data object of class `Spatial` or `sf` that contains the response variable, covariates, and coordinates
#' @param model variogram model
#' @param errorVar variances of measurement errors
#' @param min.nug minimum nugget value to avoid singular covariance matrix
#'
#' @references Hsu CC, Tsai DR, Su SY, Jhuang JR, Chiang CJ, Yang YW, Lee WC.
#'     A Stabilized Kriging Method for Mapping Disease Rates.
#'     Journal of Epidemiology (2023); 33(4):201-208.
#'
#' @export
#'
krige_sp <- function(formula, data, model, errorVar, min.nug = 0) {

  cal <- match.call()
  ind <- match(c("formula", "data", "errorVar"), names(cal), 0L)
  mf <- cal[c(1L, ind)]
  mf$na.action <- "na.fail"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  if(!inherits(data, "Spatial")) data <- as(data, "Spatial")
  data <- as(data, sub("DataFrame$", "", class(data)))
  dmat <- stats::model.matrix(formula, mf)
  nug <- pmax(get_nugget(model), min.nug)
  n <- length(data)

  d <- sp::spDists(data)
  d[d == 0] <- near0
  covar <- gstat::variogramLine(model, dist_vector = d, covariance = TRUE)
  var <- stats::model.extract(mf, "errorVar")
  if(is.null(var)) var <- nug

  C <- rbind(cbind(covar + diag(var, n), dmat),
             cbind(t(dmat), diag(0, ncol(dmat))))
  invC <- matInv(C)

  res <- list(call = cal, model = model, data = data,
              y = stats::model.response(mf), invCov = invC)
  class(res) <- c("spkrige", class(res))
  return(res)
}

#' Prediction for Spatial Kriging
#'
#' @param object fitted spatial Kriging model
#' @param newdata object of class `Spatial` or `sf` with prediction locations; should contain attributes with the covariates (if any)
#' @param computeVar logical; compute prediction variances or not
#' @param ... arguments passed to or from other methods
#'
#' @method predict spkrige
#' @export predict.spkrige
#' @export
#'
predict.spkrige <- function(object, newdata, computeVar = FALSE, ...) {

  formula <- as.formula(object$call$formula)
  data <- object$data
  model <- object$model
  n <- length(data)

  if(!inherits(newdata, "Spatial")) newdata <- as(newdata, "Spatial")
  dmat.new <- stats::model.matrix(stats::update(formula, NULL ~ .), newdata)
  d0 <- sp::spDists(data, newdata)
  d0[d0 == 0] <- near0
  covar0 <- gstat::variogramLine(model, dist_vector = d0, covariance = TRUE)
  D <- rbind(covar0, t(dmat.new))
  w <- matMult(object$invCov, D)

  # if(!all(abs(colSums(w[1:n, , drop = FALSE]) - 1) < 1e-5))
  #   warning("The kriging weights do not sum to 1")

  res <- sp::SpatialPointsDataFrame(
    as(newdata, sub("DataFrame$", "", class(newdata))),
    data = data.frame(var1.pred = c(matMult(t(object$y), w[1:n, , drop = FALSE])),
                      row.names = row.names(newdata))
  )

  if(computeVar) {
    res$var1.var <- get_sill(model) - colSums(w * D)
  }

  return(res)
}
