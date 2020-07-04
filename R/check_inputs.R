.check_input_data_matrix <- function(x, y, na.rm) {
  if (is.null(y)) {
    if (!is.array(x) & !is.data.frame(x) & !is.matrix(x)) {
      stop("\"x\" must be a numeric matrix or any other array-type object that can be converted or a \"matrix\" object with with ncol>=2.\n\n")
    }
    if (is.vector(x)) {
      stop("\"x\" must be a numeric matrix or any other array-type object that can be converted or a \"matrix\" object.  with with ncol>=2.\n\n")
    }
    if (!is.matrix(x)) {
      x <- data.matrix(x)
    }
    if (!is.numeric(x)) {
      stop("\"x\" must be numeric.")
    }
    if (nrow(x) < 2 | ncol(x) < 2) {
      stop("nrow(x)>=2 and ncol(xa)>=2 are required\n\n")
    }
    is_na_data <- is.na(x)
    if (any(is_na_data)) {
      if (na.rm == FALSE) {
        stop("\"x\" contains NA records. You may want to filter NAs by setting \"na.rm=TRUE\" (see documentation for more details).\n\n")
      }
      else {
        idx_na <- which(rowSums(is_na_data) >= 1)
        x <- x[-idx_na, , drop = FALSE]
        if (nrow(x) < 2) {
          stop("nrow(x)<2 after NA removal.\n\n")
        }
      }
    }
    if (any(!is.finite(x))) {
      stop("\"x\" contains Inf values\n\n")
    }
  }
  else {
    if (!is.vector(x) | !is.vector(y)) {
      stop("If \"y\" is given, \"x\" and  \"y\" must be both numeric.\n\n")
    }
    if (!is.numeric(x) | !is.numeric(y)) {
      stop("\"x\" and \"y\" must be numeric.\n\n")
    }
    if (length(x) != length(y)) {
      stop("\"x\" and \"y\" have different length.\n\n")
    }
    x <- cbind(x, y, deparse.level = 0)
    is_na_data <- is.na(x)
    if (any(is_na_data)) {
      if (na.rm == FALSE) {
        stop("\"x\" or \"y\" contains NA records. You may want to filter NAs by setting \"na.rm=TRUE\" (see documentation for more details).\n\n")
      }
      else {
        idx_na <- which(rowSums(is_na_data) >= 1)
        x <- x[-idx_na, , drop = FALSE]
        if (nrow(x) < 2) {
          stop("length(x)<2 and/or length(y)<2 after NA removal.\n\n")
        }
      }
    }
    if (any(!is.finite(x))) {
      stop("\"x\" and/or \"y\"  contains Inf values\n\n")
    }
  }
  return(x)
}
