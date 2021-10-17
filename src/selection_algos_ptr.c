#include <stdio.h>
/*  Implements functions quickselect and introselect */

static double select_corner_cases(double *vector_shorter_3, int size_of_vector,
                                  int extract_this_element_Ccorrected);

static double pivot_of_3(double *vector_to_be_pivoted, int size);

static int quickselect_onepass(double *vector, int size);

/* Functions */

double quickselect_recursive(double *vector, int size, int k) {
  /* Arguments:
   * *vector: pointer to vector
   * idN: size of vector
   * k: k-th element to be extracted from the vector output: the
   * extracted element of order k */

  /* Decide if we need to go with quickselect_recursive or return output
   * (quickselect_recursive works only for vectors with at least 3 elems) */
  if (size < 3) { // select the returning value
    return select_corner_cases(vector, size, k);
  }

  int subst = quickselect_onepass(vector, size);

  /* Determine where to go next (left of subst or right?) */
  if (subst == k) {
    return *(vector + subst);
  } else if (subst > k) { // go left
    size = subst;
    return quickselect_recursive(vector, size, k);
  } else {          // go right
                    // readjust k
    k -= subst + 1; // take it back to orginal value (for full vector)
    vector += subst + 1;
    size -= subst + 1;
    return quickselect_recursive(vector, size, k);
  }
}

static double select_corner_cases(double *vector, int size, int k) {
  /* Used to handel corner cases not handeled by recursive strategy.
   * Recursive strategy needs vectors of at least 3 elements

   *  Arguments:
   * *vector: pointer to vector of doubles to sort
   * size: size of the vector
   * k: k-th element to be extracted from the vector output: the
   * extracted element of order k
  */

  // total elements in vector are idN-id0+1
  double ret = -111;
  switch (size) {
  case 1:
    return *vector;
  case 2:
    switch (k) {
    case 0: // return the smallest of the two
      if (*vector < *(vector + 1))
        return *vector;
      else
        return *(vector + 1);
    case 1: // return the biggest of the two
      if (*vector > *(vector + 1))
        return *vector;
      else
        return *(vector + 1);
    }
  }
  return ret;
}

static double pivot_of_3(double *vector, int size) {
  /* perform pivot of 3 returning the median element and
   * rearranging the vector so to have:
   * id0, i(see below) , idN -> min, max, med

   * Arguments
   * vector: vector where to make substitutions
   * id0, i, idN: positions of first, median and last element to consider
  */

  int idN = size - 1; // index of last element of the vector
  int i = idN / 2;    // index of middle element of the vector
  double a, b, c;     // auxiliary variables (avoid dereferncing too much)
  a = *vector;
  b = *(vector + i);
  c = *(vector + idN);

  /* Nota; vogliamo che l'elemento mediano si trovi alla fine */
  double swapper;
  if ((a > b) ^ (a > c)) { // id0 is median element
    swapper = c;
    c = a;
    if (swapper > b) {
      a = b;
      b = swapper;
    } else {
      a = swapper;
    }
  } else if ((b > a) ^ (b > c)) { // i is median element
    swapper = c;
    c = b;
    if (swapper > a) {
      b = swapper;
    } else {
      b = a;
      a = swapper;
    }
  } else { // idN is median element
    if (a > b) {
      swapper = a;
      a = b;
      b = swapper;
    }
  }
  *vector = a;
  *(vector + i) = b;
  *(vector + idN) = c;
  return c;
}

static int quickselect_onepass(double *vector, int size) {
  /* perform one pass for quickselect and returns the element where the pivot
   * was substituted: this is used to decide whether to go left or right */

  /* Pivoting section (pivot of 3) */
  double pivot = pivot_of_3(vector, size);

  /* One pass on vector */
  double swapper;
  double *subst_el = vector; // address of element to substitute (copy not to
                             // alter main pointer)
  int subst = 0;
  for (int i = 0; i < size - 1; i++) {
    if (*subst_el < pivot) {
      if (subst == i) {
        subst++;
      } else {
        swapper = *subst_el;
        *subst_el = *(vector + subst);
        *(vector + subst) = swapper;
        subst++;
      }
    }
    subst_el++;
  }
  *subst_el = *(vector + subst);
  *(vector + subst) = pivot;
  return subst;
}
