#' Design Matrix of Age-Period-Cohort Model
#'
#' Generate the APC design matrix using the \code{designmatrix()} function from \href{https://analysistools.cancer.gov/apc/}{APC Analysis Web Tool}.
#'     The source code can be found in its \href{https://github.com/CBIIT/nci-webtools-dceg-age-period-cohort}{Github repository}.
#'     The quadratic design matrix is implemented by Chernyavskiy et al. (2019).
#'
#' @param A Age
#' @param P Period
#' @param C Cohort
#' @param degree 1 for linear APC model; 2 for quadratic APC model
#'
#' @return A matrix object
#'
#' @references Philip S. Rosenberg, David P. Check, William F. Anderson.
#'     A Web Tool for Age–Period–Cohort Analysis of Cancer Incidence and Mortality Rates.
#'     Cancer Epidemiology, Biomarkers & Prevention (2014); 23(11):2296–2302.
#' @references Chernyavskiy P, Little MP, Rosenberg PS.
#'     A unified approach for assessing heterogeneity in age–period–cohort model parameters using random effects.
#'     Statistical Methods in Medical Research (2019); 28(1):20-34.
#'
#' @export
#'
apc_dmat <- function(A, P, C, degree = 1) {
  XMAT <- cbind(A, P, C)
  PVP <- list(D = list(DATA = XMAT,
                       a = unique(XMAT[, 1]),
                       p = unique(XMAT[, 2]),
                       c = unique(XMAT[, 3])),
              RVals = apply(XMAT, 2, mean))

  if(degree == 1) {
    # A-2 age devs, P-2 period devs, C-2 cohort devs
    mat.list <- apc_dmat_linear(PVP)
  } else if(degree == 2) {
    # A-3 age devs, P-3 period devs, C-3 cohort devs
    mat.list <- apc_dmat_poly(PVP)
  }

  return(mat.list$X)
}

apc_dmat_linear <- function(PVP) {

  N <- nrow(PVP$D$DATA)
  J <- matrix(1, nrow = N)

  a <- matrix(PVP$D$DATA[, 1])
  avals <- PVP$D$a
  aM <- length(avals)


  p <- matrix(PVP$D$DATA[, 2])
  pvals <- PVP$D$p
  pN <- length(pvals)

  c <- matrix(PVP$D$DATA[, 3])
  cvals <- PVP$D$c
  cK <- length(cvals)

  # Age
  Ad <- matrix(NaN, nrow = N, ncol = aM)
  for (i in 1:aM) {
    Ad[, i] <- a == avals[i]
  }

  # Period
  Pd <- matrix(NaN, nrow = N, ncol = pN)
  for (i in 1:pN) {
    Pd[, i] <- p == pvals[i]
  }

  # Cohort
  Cd <- matrix(NaN, nrow = N, ncol = cK)
  for (i in 1:cK) {
    Cd[, i] <- c == cvals[i]
  }

  abar <- PVP$RVals[1]
  pbar <- PVP$RVals[2]
  cbar <- PVP$RVals[3]

  a0 <- a - abar
  p0 <- p - pbar
  c0 <- c - cbar

  Xa <- cbind(J, a0)
  Ra <- solve(t(Xa) %*% Xa, t(Xa))
  XAD <- diag(N) - Xa %*% Ra

  Xp <- cbind(J, p0)
  Rp <- solve(t(Xp) %*% Xp, t(Xp))
  XPD <- diag(N) - Xp %*% Rp

  Xc <- cbind(J, c0)
  Rc <- solve(t(Xc) %*% Xc, t(Xc))
  XCD <- diag(N) - Xc %*% Rc

  Ad0 <- XAD %*% Ad
  Pd0 <- XPD %*% Pd
  Cd0 <- XCD %*% Cd

  X <- cbind(J, a0, c0, Ad0[, 2:(aM-1)], Pd0[, 2:(pN-1)], Cd0[, 2:(cK-1)])

  pA <- ncol(Ad)
  pP <- ncol(Pd)
  pC <- ncol(Cd)

  Pt <- vector("list", 6)
  Pt[[1]] <- 1
  Pt[[2]] <- 2
  Pt[[3]] <- 3
  Pt[[4]] <- 4:(pA+1)
  Pt[[5]] <- (pA+2):(pA+pP-1)
  Pt[[6]] <- (pA+pP):(pA+pP+pC-3)

  # Compute contrast matrices to convert parameters to deviations
  Xa <- cbind(matrix(1, nrow = aM), matrix(avals) - abar)
  Ra <- solve(t(Xa) %*% Xa, t(Xa))
  XAD <- diag(aM) - Xa %*% Ra
  XAD <- matrix(XAD[, 2:(aM-1)], nrow = aM)

  Xp <- cbind(matrix(1, nrow = pN), matrix(pvals) - pbar)
  Rp <- solve(t(Xp) %*% Xp, t(Xp))
  XPD <- diag(pN) - Xp %*% Rp
  XPD <- matrix(XPD[, 2:(pN-1)], nrow = pN)

  Xc <- cbind(matrix(1, nrow = cK), matrix(cvals) - cbar)
  Rc <- solve(t(Xc) %*% Xc, t(Xc))
  XCD <- diag(cK) - Xc %*% Rc
  XCD <- matrix(XCD[, 2:(cK-1)], nrow = cK)
  D <- list(X = X, Pt = Pt, XAD = XAD, XPD = XPD, XCD = XCD)
  return(D)
}


apc_dmat_poly <- function(PVP) {

  N <- nrow(PVP$D$DATA)
  J <- matrix(1, nrow = N)

  a <- matrix(PVP$D$DATA[,1])
  avals <- PVP$D$a
  aM <- length(avals)

  p <- matrix(PVP$D$DATA[,2])
  pvals <- PVP$D$p
  pN <- length(pvals)

  c <- matrix(PVP$D$DATA[,3])
  cvals <- PVP$D$c
  cK <- length(cvals)

  # Age
  Ad <- matrix(NaN, nrow = N, ncol = aM)
  for (i in 1:aM) {
    Ad[, i] <- a == avals[i]
  }

  # Period
  Pd <- matrix(NaN, nrow = N, ncol = pN)
  for (i in 1:pN) {
    Pd[, i] <- p == pvals[i]
  }

  # Cohort
  Cd <- matrix(NaN, nrow = N, ncol = cK)
  for (i in 1:cK) {
    Cd[, i] <- c == cvals[i]
  }

  abar <- PVP$RVals[1]
  pbar <- PVP$RVals[2]
  cbar <- PVP$RVals[3]

  a0 <- a - abar
  p0 <- p - pbar
  c0 <- c - cbar

  a1 <- a0^2
  p1 <- p0^2
  c1 <- c0^2

  Xa <- cbind(J, a0, a1)
  Ra <- solve(t(Xa) %*% Xa, t(Xa))
  XAD <- diag(N) - Xa %*% Ra

  Xp <- cbind(J, p0, p1)
  Rp <- solve(t(Xp) %*% Xp, t(Xp))
  XPD <- diag(N) - Xp %*% Rp

  Xc <- cbind(J, c0, c1)
  Rc <- solve(t(Xc) %*% Xc, t(Xc))
  XCD <- diag(N) - Xc %*% Rc

  Ad0 <- XAD %*% Ad
  Pd0 <- XPD %*% Pd
  Cd0 <- XCD %*% Cd

  X <- cbind(J, a0, c0, a1, p1, c1, Ad0[,3:(aM-1)], Pd0[,3:(pN-1)], Cd0[,3:(cK-1)])

  pA <- ncol(Ad)
  pP <- ncol(Pd)
  pC <- ncol(Cd)

  Pt <- vector("list", 9)
  Pt[[1]] <- 1
  Pt[[2]] <- 2
  Pt[[3]] <- 3
  Pt[[4]] <- 4
  Pt[[5]] <- 5
  Pt[[6]] <- 6
  Pt[[7]] <- 7:(pA+3)
  Pt[[8]] <- (pA+4):(pA+pP)
  Pt[[9]] <- (pA+pP+1):(pA+pP+pC-3)

  # Compute contrast matrices to convert parameters to deviations
  Xa <- cbind(matrix(1, nrow = aM), matrix(avals) - abar, matrix(avals^2) - mean(avals^2))
  Ra <- solve(t(Xa) %*% Xa, t(Xa))
  XAD <- diag(aM) - Xa %*% Ra
  XAD <- matrix(XAD[, 3:(aM-1)], nrow = aM)

  Xp <- cbind(matrix(1, nrow = pN), matrix(pvals) - pbar,matrix(pvals^2) - mean(pvals^2))
  Rp <- solve(t(Xp) %*% Xp, t(Xp))
  XPD <- diag(pN) - Xp %*% Rp
  XPD <- matrix(XPD[, 3:(pN-1)], nrow = pN)

  Xc <- cbind(matrix(1, nrow = cK), matrix(cvals) - cbar, matrix(cvals^2) - mean(cvals^2))
  Rc <- solve(t(Xc) %*% Xc, t(Xc))
  XCD <- diag(cK) - Xc %*% Rc
  XCD <- matrix(XCD[, 3:(cK-1)], nrow = cK)
  D <- list(X = X, Pt = Pt, XAD = XAD, XPD = XPD, XCD = XCD)
  return(D)
}
