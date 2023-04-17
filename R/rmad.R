rmad <- function(x, y = NULL, na.rm = FALSE, even.correction = FALSE, num.threads = "half-max") {

  ## check input data
  dat <- .check_input_data_matrix(x = x, y = y, na.rm = na.rm)
  colnames_original <- colnames(dat)
  storage.mode(dat) <- "double"
  n <- as.integer(nrow(dat))
  p <- as.integer(ncol(dat))


  ## set even correction
  if (even.correction) {
    evencorrection <- 1L
  } else {
    evencorrection <- 0L
  }

  ## set number of threads
  if (num.threads == "half-max") {
    num.threads <- 0L
  } else if (num.threads == 0) {
    num.threads <- 1L
  } else {
    storage.mode(num.threads) <- "integer"
  }


  ## Call C code
  u <- .Call(C_cormad_C, dat, n, p, evencorrection, num.threads)


  if (!is.null(y)) { ## 2-dimensional
    return(u)
  } else { ## p-dimensional

    ## assemble the matrix using the lower triangle
    R <- Matrix(1, nrow = p, ncol = p, sparse = FALSE)
    R[lower.tri(R, diag = FALSE)] <- u
    R <- forceSymmetric(R, uplo = "L")
    R <- as(R, "packedMatrix")


    ## attach dimnames if needed
    if (!is.null(colnames_original)) {
      dimnames(R)[[1]] <- dimnames(R)[[2]] <- colnames_original
    }

    return(R)
  } ## END if(!is.null(y)){ ## 2-dimensional
} ## END function
