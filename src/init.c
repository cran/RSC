#include <stdlib.h>
//#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/RS.h>
#include <R_ext/Rdynload.h>

extern void F77_NAME(cormadvecdp)(void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortMethods[] = {
   {"cormadvecdp", (DL_FUNC) &F77_NAME(cormadvecdp), 6},
   {NULL, NULL, 0}
};

void
R_init_RSC(DllInfo *dll)
{
   R_registerRoutines(dll, NULL, NULL, FortMethods, NULL);
   R_useDynamicSymbols(dll, FALSE);
}
