#include "RSCdefines.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

// Wrapper for C cormad (with pointers)
SEXP cormad_C(SEXP R_matrix, SEXP len_rows, SEXP len_cols, SEXP R_evencorrect,
              SEXP R_num_threads) {
  // prepare input arguments
  int n_row = asInteger(len_rows);
  int n_col = asInteger(len_cols);
  int output_size = (n_col - 1) * n_col / 2;
  int evencorrect =
      (n_row % 2 == 1) ? 0 : asInteger(R_evencorrect); // ingore for odd n
  int num_threads = asInteger(R_num_threads);
  SEXP output = PROTECT(allocVector(REALSXP, output_size)); // define output

  // create pointers to (duplicate) matrix column
  double *dupmat =
      REAL(PROTECT(duplicate(R_matrix))); // copy input (to avoid modifications)

  // Call cormad
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (num_threads == 0) /*  use max number of threads / 2 */
    num_threads = (max_threads / 2 == 0) ? 1 : max_threads / 2;
  else if (num_threads > max_threads)
    num_threads = max_threads;
  else if (num_threads < 0)
    num_threads =
        (max_threads + num_threads > 0) ? max_threads + num_threads : 1;

  if (num_threads > 1) {
    cormad_parallel(dupmat, n_row, n_col, REAL(output), evencorrect,
                    num_threads);
  } else {
    cormad(dupmat, n_row, n_col, REAL(output), evencorrect);
  }
#else
  if (num_threads != 1)
    Rprintf("\nSpecified num_threads=%d, but OPENMP support not found; "
            "switching to single core.\n\n",
            num_threads);
  cormad(dupmat, n_row, n_col, REAL(output), evencorrect);
#endif
  // Remove protect to allow R grabage collect
  UNPROTECT(2);
  return output;
}

static const R_CallMethodDef callMethods[] = {
    {"cormad_C", (DL_FUNC)&cormad_C, 5}, {NULL, NULL, 0}};

void R_init_RSC(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
