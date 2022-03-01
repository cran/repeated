#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* .C calls */
extern void countfb_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ddb_c(void *, void *, void *, void *, void *, void *, void *);
extern void ddp_c(void *, void *, void *, void *, void *, void *, void *);
extern void dmb_c(void *, void *, void *, void *, void *, void *, void *);
extern void dmp_c(void *, void *, void *, void *, void *, void *, void *);
extern void dpvfp_c(void *, void *, void *, void *, void *, void *, void *);
extern void gar_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kcountb_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void krand_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kserie_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pginvgauss_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ppowexp_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void psimplex_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
//extern void romberg_c(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"countfb_c",    (DL_FUNC) &countfb_c,    21},
    {"ddb_c",        (DL_FUNC) &ddb_c,         7},
    {"ddp_c",        (DL_FUNC) &ddp_c,         7},
    {"dmb_c",        (DL_FUNC) &dmb_c,         7},
    {"dmp_c",        (DL_FUNC) &dmp_c,         7},
    {"dpvfp_c",      (DL_FUNC) &dpvfp_c,       7},
    {"gar_c",        (DL_FUNC) &gar_c,        25},
    {"kcountb_c",    (DL_FUNC) &kcountb_c,    24},
    {"krand_c",      (DL_FUNC) &krand_c,      25},
    {"kserie_c",     (DL_FUNC) &kserie_c,     26},
    {"pginvgauss_c", (DL_FUNC) &pginvgauss_c, 10},
    {"ppowexp_c",    (DL_FUNC) &ppowexp_c,    10},
    {"psimplex_c",   (DL_FUNC) &psimplex_c,   10},
//   {"romberg_c",    (DL_FUNC) &romberg_c,     9},
    {NULL, NULL, 0}
};


/* .Call()  */
//extern void romberg_sexp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
SEXP romberg_sexp(SEXP fcn, SEXP a, SEXP b, SEXP len, SEXP eps,
                   SEXP pts, SEXP max, SEXP err, SEXP envir);

static const R_CallMethodDef callMethods[]  = {
  {"romberg_sexp", (DL_FUNC) &romberg_sexp, 9},
  {NULL, NULL, 0}
};


/* .Fortran calls */
extern void F77_NAME(binnest_f)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(chidden_f)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cphidden_f)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gcopula_f)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hidden_f)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logitord_f)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"binnest_f",  (DL_FUNC) &F77_NAME(binnest_f),  64},
    {"chidden_f",  (DL_FUNC) &F77_NAME(chidden_f),  40},
    {"cphidden_f", (DL_FUNC) &F77_NAME(cphidden_f), 40},
    {"gcopula_f",  (DL_FUNC) &F77_NAME(gcopula_f),  17},
    {"hidden_f",   (DL_FUNC) &F77_NAME(hidden_f),   35},
    {"logitord_f", (DL_FUNC) &F77_NAME(logitord_f), 19},
    {NULL, NULL, 0}
};

void R_init_repeated(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, callMethods, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
