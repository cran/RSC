/* Defines quickselect and cormad (double precision) */

double quickselect_recursive(double *vector_extract_k, int vec_size, int extract_this_element);

void cormad(double *matrix, int n_row, int n_col, double *output, int evencorrect);

#ifdef _OPENMP
void cormad_parallel(double *matrix, int n_row, int n_col, double *output,
                      int evencorrect, int num_threads);
#endif

