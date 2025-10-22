#' Design Matrix of Age-Period-Cohort Model
#'
#' Construct the APC design matrix with an age-drift parameterization.
#'
#' @param A Age
#' @param P Period
#' @param C Cohort
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

  Apos <- outer(A, a, `==`)
  Ppos <- outer(P, p, `==`)
  Cpos <- outer(C, c, `==`)

  A0 <- A - M[1]
  P0 <- P - M[2]
  C0 <- C - M[3]

  if(degree == 2) {
    A0 <- cbind(A0, A0^2)
    P0 <- cbind(P0, P0^2)
    C0 <- cbind(C0, C0^2)
  }

  Ax <- cbind(1, A0)
  Acoef <- solve(crossprod(Ax), t(Ax))
  Adm <- (I - Ax %*% Acoef) %*% Apos

  Px <- cbind(1, P0)
  Pcoef <- solve(crossprod(Px), t(Px))
  Pdm <- (I - Px %*% Pcoef) %*% Ppos

  Cx <- cbind(1, C0)
  Ccoef <- solve(crossprod(Cx), t(Cx))
  Cdm <- (I - Cx %*% Ccoef) %*% Cpos

  if(degree == 1) {

    # A-2 age devs, P-2 period devs, C-2 cohort devs
    dm <- cbind(1, A0, C0,
                Adm[, -c(1, ncol(Adm))],
                Pdm[, -c(1, ncol(Pdm))],
                Cdm[, -c(1, ncol(Cdm))], deparse.level = 0)

  } else if(degree == 2) {

    # A-3 age devs, P-3 period devs, C-3 cohort devs
    dm <- cbind(1, A0[, 1], C0[, 1], A0[, 2], P0[, 2], C0[, 2],
                Adm[, -c(1:2, ncol(Adm))],
                Pdm[, -c(1:2, ncol(Pdm))],
                Cdm[, -c(1:2, ncol(Cdm))], deparse.level = 0)
  }

  return(dm)
}
