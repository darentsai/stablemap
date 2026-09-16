#' Spatiotemporal Kriging (Ordinary or Universal)
#'
#' @param formula formula that defines the response variable as a linear model of covariates
#' @param data object of class `ST*DF` that contains the response variable, covariates, and coordinates
#' @param model spatio-temporal variogram model
#' @param errorVar variances of measurement errors
#' @param min.nug minimum nugget value to avoid singular covariance matrix
#'
#' @references Tsai DR, Jhuang JR, Su SY, Chiang CJ, Yang YW, Lee WC.
#'     A stabilized spatiotemporal kriging method for disease mapping and application to male oral cancer and female breast cancer in Taiwan.
#'     BMC Medical Research Methodology (2022); 22(1):270.
#'
#' @export
#'
krige_st <- function(formula, data, model, errorVar, min.nug = 0) {

  cal <- match.call()
  ind <- match(c("formula", "data", "errorVar"), names(cal), 0L)
  mf <- cal[c(1L, ind)]
  mf$na.action <- "na.fail"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  data <- as(data, "STS")
  dmat <- stats::model.matrix(formula, mf)
  nug <- pmax(get_nugget(model), min.nug)
  sill <- get_sill(model)
  n <- length(data)

  ind <- data@index
  d.s <- sp::spDists(data@sp)
  d.s[d.s == 0] <- near0
  d.t <- as.matrix(dist(spacetime::index(data@time)))
  d.t[d.t == 0] <- near0
  d <- data.frame(spacelag = d.s[t(combn(ind[, 1], 2))],
                  timelag = d.t[t(combn(ind[, 2], 2))])

  covar <- array(dim = c(n, n))
  covar[lower.tri(covar)] <- gstat::variogramSurface(model, dist_grid = d, covariance = TRUE)$gamma
  covar[upper.tri(covar)] <- t(covar)[upper.tri(covar)]
  diag(covar) <- sill
  var <- stats::model.extract(mf, "errorVar")
  if(is.null(var)) var <- nug

  C <- rbind(cbind(covar + diag(var, n), dmat),
             cbind(t(dmat), diag(0, ncol(dmat))))
  invC <- matInv(C)

  res <- list(call = cal, model = model, data = data, y = stats::model.response(mf), invCov = invC)
  class(res) <- c("stkrige", class(res))
  return(res)
}

#' Prediction for Spatiotemporal Kriging
#'
#' @param object fitted spatiotemporal Kriging model
#' @param newdata object of class `ST*` with prediction locations; should contain attributes with the covariates (if any)
#' @param computeVar logical; compute prediction variances or not
#' @param ... arguments passed to or from other methods
#'
#' @method predict stkrige
#' @export predict.stkrige
#' @export
#'
predict.stkrige <- function(object, newdata, computeVar = FALSE, ...) {

  formula <- as.formula(object$call$formula)
  data <- object$data
  model <- object$model
  ind <- data@index
  n <- length(data)

  newdata <- as(newdata, sub("ST.", "STS", class(newdata)))
  dmat.new <- stats::model.matrix(stats::update(formula, NULL ~ .), newdata)
  ind0 <- newdata@index
  d0.s <- sp::spDists(data@sp, newdata@sp)
  d0.s[d0.s == 0] <- near0
  d0.t <- abs(outer(spacetime::index(data@time), spacetime::index(newdata@time), `-`))
  d0.t[d0.t == 0] <- near0
  if(is(d0.t, "difftime")) d0.t <- unclass(d0.t)
  d0 <- data.frame(spacelag = d0.s[as.matrix(expand.grid(data = ind[, 1], newdata = ind0[, 1]))],
                   timelag = d0.t[as.matrix(expand.grid(data = ind[, 2], newdata = ind0[, 2]))])
  covar0 <- matrix(gstat::variogramSurface(model, dist_grid = d0, covariance = TRUE)$gamma, nrow = n)
  D <- rbind(covar0, t(dmat.new))
  w <- matMult(object$invCov, D)

  # if(!all(abs(colSums(w[1:n, , drop = FALSE]) - 1) < 1e-5))
  #   warning("The kriging weights do not sum to 1")

  res <- spacetime::STSDF(
    sp = newdata@sp, time = newdata@time,
    index = newdata@index, endTime = newdata@endTime,
    data = data.frame(var1.pred = c(matMult(t(object$y), w[1:n, , drop = FALSE])))
  )

  if(computeVar) {
    res$var1.var <- get_sill(model) - colSums(w * D)
  }

  return(res)
}
