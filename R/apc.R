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
  I <- diag(N)

  a <- sort(unique(A))
  p <- sort(unique(P))
  c <- sort(unique(C))

  na <- length(a)
  np <- length(p)
  nc <- length(c)

  ma <- floor((na + 1) / 2)
  mp <- floor((np + 1) / 2)
  mc <- mp - ma + na

  aref <- a[ma]
  pref <- p[mp]
  cref <- c[mc]

  A0 <- matrix(A - aref)
  P0 <- matrix(P - pref)
  C0 <- matrix(C - cref)

  if(degree == 2) {
    A0 <- cbind(A0, A0^2)
    P0 <- cbind(P0, P0^2)
    C0 <- cbind(C0, C0^2)
  }

  # one-hot dummy matrix
  Ad <- outer(A, a, `==`)
  Pd <- outer(P, p, `==`)
  Cd <- outer(C, c, `==`)

  # age curvatures
  Ax <- cbind(1, A0)
  Ay <- solve(crossprod(Ax), t(Ax))
  Acurv <- (I - Ax %*% Ay) %*% Ad

  # period curvatures
  Px <- cbind(1, P0)
  Py <- solve(crossprod(Px), t(Px))
  Pcurv <- (I - Px %*% Py) %*% Pd

  # cohort curvatures
  Cx <- cbind(1, C0)
  Cy <- solve(crossprod(Cx), t(Cx))
  Ccurv <- (I - Cx %*% Cy) %*% Cd

  if(degree == 1) {
    # A-2 age curvatures, P-2 period curvatures, C-2 cohort curvatures
    dmat <- cbind(1, A0, C0,
                  Acurv[, 2:(na-1)], Pcurv[, 2:(np-1)], Ccurv[, 2:(nc-1)],
                  deparse.level = 0)

  } else if(degree == 2) {
    # A-3 age curvatures, P-3 period curvatures, C-3 cohort curvatures
    dmat <- cbind(1, A0[, 1], C0[, 1], A0[, 2], P0[, 2], C0[, 2],
                  Acurv[, 3:(na-1)], Pcurv[, 3:(np-1)], Ccurv[, 3:(nc-1)],
                  deparse.level = 0)
  }

  return(dmat)
}
