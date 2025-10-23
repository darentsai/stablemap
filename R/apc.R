#' Design Matrix of Age-Period-Cohort Model
#'
#' Construct a linear or quadratic APC design matrix with an age-drift parameterization.
#'
#' @param A Age index
#' @param P Period index
#' @param C Cohort index
#' @param degree 1 for linear APC model; 2 for quadratic APC model
#'
#' @return A matrix object
#'
#' @export
#'
apc_dmat <- function(A, P, C, degree = 1) {

  dat <- cbind(A, P, C)
  N <- nrow(dat)
  M <- colMeans(dat)
  I <- diag(N)

  a <- unique(A)
  p <- unique(P)
  c <- unique(C)

  # one-hot dummy matrix
  Adum <- outer(A, a, `==`)
  Pdum <- outer(P, p, `==`)
  Cdum <- outer(C, c, `==`)

  A0 <- A - M[1]
  P0 <- P - M[2]
  C0 <- C - M[3]

  if(degree == 2) {
    A0 <- cbind(A0, A0^2)
    P0 <- cbind(P0, P0^2)
    C0 <- cbind(C0, C0^2)
  }

  # age curvatures
  Ax <- cbind(1, A0)
  Acoef <- solve(crossprod(Ax), t(Ax))
  Acurv <- (I - Ax %*% Acoef) %*% Adum

  # period curvatures
  Px <- cbind(1, P0)
  Pcoef <- solve(crossprod(Px), t(Px))
  Pcurv <- (I - Px %*% Pcoef) %*% Pdum

  # cohort curvatures
  Cx <- cbind(1, C0)
  Ccoef <- solve(crossprod(Cx), t(Cx))
  Ccurv <- (I - Cx %*% Ccoef) %*% Cdum

  if(degree == 1) {

    # A-2 age curvs, P-2 period curvs, C-2 cohort curvs
    dmat <- cbind(1, A0, C0,
                  Acurv[, -c(1, ncol(Acurv))],
                  Pcurv[, -c(1, ncol(Pcurv))],
                  Ccurv[, -c(1, ncol(Ccurv))], deparse.level = 0)

  } else if(degree == 2) {

    # A-3 age curvs, P-3 period curvs, C-3 cohort curvs
    dmat <- cbind(1, A0[, 1], C0[, 1], A0[, 2], P0[, 2], C0[, 2],
                  Acurv[, -c(1:2, ncol(Acurv))],
                  Pcurv[, -c(1:2, ncol(Pcurv))],
                  Ccurv[, -c(1:2, ncol(Ccurv))], deparse.level = 0)
  }

  return(dmat)
}
