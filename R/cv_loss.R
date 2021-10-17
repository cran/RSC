## Normalized squared Frobenius loss for all threshold values at a given train/test
## set, where
##
##   * idx = TRUE for sample points into the train set
##   * train set = dat[  idx , ]
##   * test  set = dat[ !idx , ]
##
.cv_loss <- function(idx, dat, evencorrection, threshold, grid.length, p) {


  ## compute RMAD  on train set
  n1 <- as.integer(sum(idx))
  C1 <- .Call(C_cormad_C, dat[idx, ], n1, p, evencorrection, num.threads = 1)

  ## compute RMAD  on test set
  n2 <- as.integer(sum(!idx))
  C2 <- .Call(C_cormad_C, dat[!idx, ], n2, p, correcteven = evencorrection, num.threads = 1)

  ## apply thresholds
  ans <- rep(0, times = grid.length)
  for (h in 1:grid.length) {
    C1[abs(C1) < threshold[h]] <- 0 ## fit on train set
    ans[h] <- sum(2 * {
      C1 - C2
    }^2) / p
  }

  return(ans)
}
