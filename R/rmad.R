rmad <- function(x, y = NULL, na.rm = FALSE, even.correction = FALSE) {
  dat <- .check_input_data_matrix(x = x, y = y, na.rm = na.rm)
  colnames_original <- colnames(dat)
  storage.mode(dat) <- "double"
  n <- as.integer(nrow(dat))
  p <- as.integer(ncol(dat))
  nc <- as.integer({
    p^2 - p
  }/2)
  if (even.correction) {
    evencorrection <- 1L
  }
  else {
    evencorrection <- 0L
  }
  u <- .Fortran("cormadvecdp", matrix = dat, nrow = n, ncol = p, res = numeric(nc), 
    ressize = nc, correcteven = evencorrection, PACKAGE = "RSC")$res
  if (!is.null(y)) {
    return(u)
  }
  else {
    R <- Matrix(1, nrow = p, ncol = p, sparse = FALSE)
    R[lower.tri(R, diag = FALSE)] <- u
    R <- forceSymmetric(R, uplo = "L")
    R <- as(R, "dspMatrix")
    if (!is.null(colnames_original)) {
      dimnames(R)[[1]] <- dimnames(R)[[2]] <- colnames_original
    }
    return(R)
  }
}
