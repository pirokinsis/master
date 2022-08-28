
#include <stdlib.h> // for NULL


/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(dgelyp)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gclmll)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gclmls)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"dgelyp",    (DL_FUNC) &F77_NAME(dgelyp),     8},
    {"gclmll",    (DL_FUNC) &F77_NAME(gclmll),    11},
    {"gclmls",    (DL_FUNC) &F77_NAME(gclmls),    11},
    {NULL, NULL, 0}
};

void R_init_clggm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
