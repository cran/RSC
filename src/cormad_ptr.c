#include "RSCdefines.h"
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define SQRT2 (sqrt(2))
#define CONST 1.4826

void cormad(double *matrix, int n_row, int n_col, double *output,
            int evencorrect) {

  int k = n_row / 2; // position of the median
  int output_size = (n_col - 1) * n_col / 2;
  double med, mad;           // store variables for medians
  double U[n_row], V[n_row]; // help vectors

  /* Transform matrix columns (CORMAD part 1)
   * each column has n entries;
   * matrix is assumed to be streamed in onevector
   */
  int i = 0;
  for (int l = 0; l < n_col * n_row; l++) {
    if (i == (n_row - 1)) {
      U[i] = *matrix;
      med = quickselect_recursive(U, n_row, k);
      if (evencorrect == 1) { // handle even correction
        med = (med + quickselect_recursive(U, n_row, k - 1)) / 2;
      }
      for (int j = 0; j < n_row; j++) {
        U[j] = *(matrix - (n_row - 1) + j) - med;
        V[j] = fabs(U[j]);
      }
      med = quickselect_recursive(V, n_row, k);
      if (evencorrect == 1) { // handle even correction
        med = (med + quickselect_recursive(V, n_row, k - 1)) / 2;
      }
      for (int j = 0; j < n_row; j++) // reassign
        *(matrix - (n_row - 1) + j) = U[j] / (SQRT2 * CONST * med);
      // prepare next iter
      matrix++;
      i = 0;
    } else {
      U[i] = *matrix;
      // prepare next iter
      matrix++;
      i++;
    }
  }

  matrix = matrix - (n_col * n_row); /* reset pointer */

  /* Operate on columns pairs (CORMAD part 2) */
  double *matrix_2 = matrix; // used to point at second column
  int first_col = 0;         /* Running first column */
  int second_col = 0;        /* Running second column */
  for (int l = 0; l < output_size; l++) {
    if (second_col == n_col - 1) {
      first_col++;
      second_col = first_col;
      /* set pointers to columns */
      matrix += n_row;
      matrix_2 = matrix;
    }
    second_col++;
    matrix_2 += n_row;

    for (int i = 0; i < n_row; i++) { // auxiliary vectors from matrix
      U[i] = *(matrix + i) + *(matrix_2 + i);
      V[i] = -*(matrix + i) + *(matrix_2 + i);
    }
    mad = quickselect_recursive(U, n_row, k);
    med = quickselect_recursive(V, n_row, k);
    if (evencorrect == 1) {
      mad = (mad + quickselect_recursive(U, n_row, k - 1)) / 2;
      med = (med + quickselect_recursive(V, n_row, k - 1)) / 2;
    }
    for (int i = 0; i < n_row; i++) { // reassign for new medians
      U[i] = fabs(U[i] - mad);
      V[i] = fabs(V[i] - med);
    }
    mad = quickselect_recursive(U, n_row, k);
    med = quickselect_recursive(V, n_row, k);
    if (evencorrect == 1) {
      mad = (mad + quickselect_recursive(U, n_row, k - 1)) / 2;
      med = (med + quickselect_recursive(V, n_row, k - 1)) / 2;
    }
    mad = pow(CONST * mad, 2);
    med = pow(CONST * med, 2);

    // Assign output
    *output = (mad - med) / (mad + med);
    output++;
  }
}

#ifdef _OPENMP
void cormad_parallel(double *matrix, int n_row, int n_col, double *output,
                     int evencorrect, int num_threads) {
  int k = n_row / 2; // position of the median
  int output_size = (n_col - 1) * n_col / 2;
  double med, mad;           // store variables for medians
  double U[n_row], V[n_row]; // help vectors

  /* Transform matrix columns (CORMAD part 1) */
  double *help_matrix = matrix; /*  help pointer for matrix */
#pragma omp parallel for private(med, mad, U, V, help_matrix)                  \
    num_threads(num_threads)
  for (int j = 0; j < n_col; j++) {   // iterate on cols
    help_matrix = matrix + n_row * j; /* set pointer at beg of col */
    for (int i = 0; i < n_row; i++)
      U[i] = *(help_matrix + i);
    med = quickselect_recursive(U, n_row, k);
    if (evencorrect == 1) { // handle even correction
      med = (med + quickselect_recursive(U, n_row, k - 1)) / 2;
    }
    for (int i = 0; i < n_row; i++) {
      U[i] = *(help_matrix + i) - med;
      V[i] = fabs(U[i]);
    }
    med = quickselect_recursive(V, n_row, k);
    if (evencorrect == 1) { // handle even correction
      med = (med + quickselect_recursive(V, n_row, k - 1)) / 2;
    }
    for (int i = 0; i < n_row; i++) // reassign
      *(help_matrix + i) = U[i] / (SQRT2 * CONST * med);
  }

  // int l = 0; // used to iterate over output
  /* Operate on columns pairs (CORMAD part 2) */
  help_matrix = matrix; /* help pointers for matrix */
  double *help_matrix_2 = matrix;
#pragma omp parallel for num_threads(num_threads) private(                     \
    med, mad, U, V, help_matrix, help_matrix_2)
  for (int l = 0; l < output_size; l++) {
    /*  Detrmine columns pairs */
    int col1 = 0, col2 = 0;
    int copy_l = l + 1;
    for (int last_elem = n_col - 1; last_elem > 0; last_elem--) {
      copy_l -= last_elem;
      if (copy_l <= 0) {
        col2 = (n_col - 1) + copy_l;
        break;
      } else {
        col1++;
      }
    }
    /* Set pointers to columns */
    help_matrix = matrix + n_row * col1;   // col1
    help_matrix_2 = matrix + n_row * col2; // col2

    for (int i = 0; i < n_row; i++) { // auxiliary vectors from matrix
      U[i] = *(help_matrix + i) + *(help_matrix_2 + i);
      V[i] = -*(help_matrix + i) + *(help_matrix_2 + i);
    }
    mad = quickselect_recursive(U, n_row, k);
    med = quickselect_recursive(V, n_row, k);
    if (evencorrect == 1) {
      mad = (mad + quickselect_recursive(U, n_row, k - 1)) / 2;
      med = (med + quickselect_recursive(V, n_row, k - 1)) / 2;
    }
    for (int i = 0; i < n_row; i++) { // reassign for new medians
      U[i] = fabs(U[i] - mad);
      V[i] = fabs(V[i] - med);
    }
    mad = quickselect_recursive(U, n_row, k);
    med = quickselect_recursive(V, n_row, k);
    if (evencorrect == 1) {
      mad = (mad + quickselect_recursive(U, n_row, k - 1)) / 2;
      med = (med + quickselect_recursive(V, n_row, k - 1)) / 2;
    }
    mad = pow(CONST * mad, 2);
    med = pow(CONST * med, 2);

    // Assign output
    *(output + l) = (mad - med) / (mad + med);
  }
}
#endif
