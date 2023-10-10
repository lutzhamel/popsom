#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(vsom)(float *neurons,
                           const float *dt,
                           const int *dtrows,
                           const int *dtcols,
                           const int *xdim,
                           const int *ydim,
                           const float *alpha,
                           const int *train,
                           const int *seed);

static const R_FortranMethodDef FortranEntries[] = {
    {"vsom", (DL_FUNC) &F77_NAME(vsom), 9},
    {NULL, NULL, 0}
};

void R_init_popsom(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
