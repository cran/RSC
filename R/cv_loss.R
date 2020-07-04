.cv_loss <- function(idx, dat, evencorrection, threshold, grid.length, p, nc) {
  res <- numeric(nc)
  n1 <- as.integer(sum(idx))
  C1 <- .Fortran("cormadvecdp", matrix = dat[idx, ], nrow = n1, ncol = p, res = res, 
    ressize = nc, correcteven = evencorrection, PACKAGE = "RSC")$res
  n2 <- as.integer(sum(!idx))
  C2 <- .Fortran("cormadvecdp", matrix = dat[!idx, ], nrow = n2, ncol = p, res = res, 
    ressize = nc, correcteven = evencorrection, PACKAGE = "RSC")$res
  ans <- rep(0, times = grid.length)
  for (h in 1:grid.length) {
    C1[abs(C1) < threshold[h]] <- 0
    ans[h] <- sum(2 * {
      C1 - C2
    }^2)/p
  }
  return(ans)
}
