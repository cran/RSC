rsc_cv <- function(x, cv.type = "kfold", R = 10, K = 10, threshold = seq(0.05, 0.95, 
  by = 0.025), even.correction = FALSE, na.rm = FALSE, ncores = NULL, monitor = TRUE) {
  dat <- .check_input_data_matrix(x = x, y = NULL, na.rm = na.rm)
  colnames_original <- colnames(dat)
  storage.mode(dat) <- "double"
  n <- as.integer(nrow(dat))
  p <- as.integer(ncol(dat))
  nc <- as.integer({
    p^2 - p
  }/2)
  if ({
    cv.type != "random"
  } & {
    cv.type != "kfold"
  }) {
    stop("\"cv.type\" must be either \"random\" (default)  or \"kfold\"")
  }
  if (!is.numeric(R)) {
    stop("\"R\" must be an integer > 1")
  }
  else if (R < 1) {
    stop("\"R\" must be an integer > 1")
  }
  if (!is.numeric(K)) {
    stop("\"K\" must be an integer > 1")
  }
  else if (R < 1) {
    stop("\"K\" must be an integer > 1")
  }
  if (length(threshold) == 1) {
    if (threshold <= 0 | threshold >= 1) {
      stop("\"threshold\" value does not belong to the interval (0,1).")
    }
  }
  else if (length(threshold) > 1) {
    if (any(threshold < 0) | any(threshold >= 1)) {
      stop("Some of the \"threshold\" values do not belong to the interval (0,1).")
    }
    grid.length <- length(threshold)
  }
  if (even.correction) {
    evencorrection <- 1L
  }
  else {
    evencorrection <- 0L
  }
  if (is.null(ncores)) {
    DetectedCores <- detectCores()
    if (DetectedCores <= 2) {
      ncores <- 1
    }
    else {
      ncores <- {
        DetectedCores - 1
      }
    }
  }
  else {
    ncores <- as.integer(ncores)
    if (ncores <= 0) {
      stop("\"ncores\" must be an integer larger or equal to 1.")
    }
  }
  if (monitor) {
    cat("\n")
    message("Computing the RMAD matrix")
    t0 <- Sys.time()
  }
  rmad_vec <- .Fortran("cormadvecdp", matrix = dat, nrow = n, ncol = p, res = numeric(nc), 
    ressize = nc, correcteven = evencorrection, PACKAGE = "RSC")$res
  if (monitor) {
    t1 <- Sys.time()
    dt01 <- difftime(t1, t0, units = "auto")
    message("* RMAD computing time:...... ", round(dt01, 2), " [", attributes(dt01)$units, 
      "]")
  }
  if (cv.type == "random") {
    if (monitor) {
      cat("\n")
      message("Performing cross-validation")
      t_hat <- round(1.2 * {
        {
          dt01 * R * 2
        }/ncores
      }, 2)
      message("* predicted end time (worst case):...... ", Sys.time() + t_hat)
    }
    n1 <- n - floor(n/log(n))
    IDX <- array(FALSE, dim = c(R, n))
    for (r in 1:R) {
      IDX[r, ][sample(1:n, size = n1, replace = FALSE)] <- TRUE
    }
    registerDoParallel(ncores)
    U <- foreach(r = 1:R) %dopar% {
      .cv_loss(idx = IDX[r, ], dat = dat, evencorrection = evencorrection, 
        threshold = threshold, grid.length = grid.length, p = p, nc = nc)
    }
    stopImplicitCluster()
    LOSS <- array(0, dim = c(R, grid.length))
    for (r in 1:R) {
      LOSS[r, ] <- U[[r]]
    }
    avg_loss <- apply(LOSS, 2, mean)
    se_loss <- apply(LOSS, 2, sd)/sqrt(R)
  }
  if (cv.type == "kfold") {
    if (monitor) {
      cat("\n")
      message("Performing cross-validation")
      t_hat <- round(1.2 * {
        {
          dt01 * R * K * 2
        }/ncores
      }, 2)
      message("* predicted end time (worst case):...... ", Sys.time() + t_hat)
    }
    idx_fold <- cut(1:n, breaks = K, labels = FALSE)
    IDX <- array(TRUE, dim = c(R * K, n))
    row_count <- 1L
    for (r in 1:R) {
      idx_fold_shuffle <- sample(idx_fold, size = n, replace = FALSE)
      for (k in 1:K) {
        IDX[row_count, ][idx_fold_shuffle == k] <- FALSE
        row_count <- 1L + row_count
      }
    }
    registerDoParallel(ncores)
    U <- foreach(r = 1:{
      R * K
    }) %dopar% {
      .cv_loss(idx = IDX[r, ], dat = dat, evencorrection = evencorrection, 
        threshold = threshold, grid.length = grid.length, p = p, nc = nc)
    }
    stopImplicitCluster()
    if (R == 1) {
      LOSS <- array(0, dim = c(K, grid.length))
      for (k in 1:K) {
        LOSS[k, ] <- U[[k]]
      }
    }
    else {
      LOSS <- array(0, dim = c(K, grid.length, R))
      dimnames(LOSS)[[3]] <- paste0("r", 1:R)
      dimnames(LOSS)[[1]] <- paste0("k", 1:K)
      row_count <- 1L
      for (r in 1:R) {
        for (k in 1:K) {
          LOSS[k, , r] <- U[[row_count]]
          row_count <- 1L + row_count
        }
      }
    }
    avg_loss_r <- sd_loss_r <- matrix(0, nrow = R, ncol = grid.length)
    for (r in 1:R) {
      avg_loss_r[r, ] <- apply(LOSS[, , r], 2, mean)
      sd_loss_r[r, ] <- apply(LOSS[, , r], 2, sd)/sqrt(K)
    }
    avg_loss <- apply(avg_loss_r, 2, mean)
    se_loss <- apply(sd_loss_r, 2, mean)
    if (monitor) {
      t2 <- Sys.time()
      dt02 <- difftime(t2, t0, units = "auto")
    }
  }
  if (monitor) {
    message("* finished on:.......................... ", Sys.time())
  }
  tstar <- which.min(avg_loss)
  flags <- rep("", grid.length)
  a <- avg_loss[tstar] - se_loss[tstar]
  b <- avg_loss[tstar] + se_loss[tstar]
  flags[{
    avg_loss >= a
  } & {
    avg_loss <= b
  }] <- "*"
  flags[tstar] <- "minimum"
  res <- data.frame(Threshold = threshold, Average = avg_loss, SE = se_loss, Flag = flags)
  ans <- list(rmadvec = rmad_vec, varnames = colnames_original, loss = res, minimum = threshold[tstar], 
    minimum1se = max(threshold[avg_loss >= a & avg_loss <= b]))
  class(ans) <- "rsc_cv"
  return(ans)
}
