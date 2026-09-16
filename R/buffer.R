#' Create Circular Buffers
#'
#' @param data data for which buffers need to be added
#' @param ref data used as a reference to generate buffers
#' @param buffvar variable for which buffers are created
#' @param timevar time variable
#' @param sdist,tdist buffer sizes for space (in meters) and time
#'
#' @return the design matrix for buffers
#'
#' @export
#'
buffer <- function(data, ref, buffvar, timevar, sdist, tdist) {

  dat0 <- data
  if(!is(data, "sf")) data <- sf::st_as_sf(data)
  if(!is(ref, "sf")) ref <- sf::st_as_sf(ref)

  g1 <- sf::st_geometry(data); g1x <- unique(g1)
  g2 <- sf::st_geometry(ref); g2x <- unique(g2)
  dmat <- sf::st_distance(g1x, g2x)
  t1 <- data[[timevar]]; t1x <- unique(t1)
  t2 <- ref[[timevar]]
  indx <- cbind(match(g1, g1x), match(t1, t1x))

  res <- mapply(\(ds, dt) {
    buff <- sapply(t1x, \(t) {
      refx <- ref[t2 >= t - dt & t2 <= t + dt, ]
      pos <- match(st_geometry(refx), g2x)
      apply(dmat, 1, \(d) {
        mean(refx[[buffvar]][d[pos] < ds], na.rm = TRUE)
      })
    })
    return(buff[indx])
  }, sdist, tdist)

  buffname <- paste0(buffvar, "_buffer", seq_len(ncol(res)))
  dat0[buffname] <- res[ , , drop = TRUE]
  return(dat0)
}
