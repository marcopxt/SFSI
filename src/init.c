/*
 This file was autogenerated with the following R routine:
 tools::package_native_routine_registration_skeleton("package_SFSI/SFSI")
*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cov2correlation(SEXP, SEXP, SEXP, SEXP);
extern SEXP cov2distance(SEXP, SEXP, SEXP);
extern SEXP delete_col(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getCorrelated(SEXP, SEXP, SEXP);
extern SEXP readBinFileFloat(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP updatebeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP writeBinFileFloat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP addvalue2diag(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cov2correlation",   (DL_FUNC) &cov2correlation,    4},
    {"cov2distance",      (DL_FUNC) &cov2distance,       3},
    {"delete_col",        (DL_FUNC) &delete_col,         6},
    {"getCorrelated",     (DL_FUNC) &getCorrelated,      3},
    {"readBinFileFloat",  (DL_FUNC) &readBinFileFloat,   5},
    {"updatebeta",        (DL_FUNC) &updatebeta,        11},
    {"writeBinFileFloat", (DL_FUNC) &writeBinFileFloat,  6},
    {"addvalue2diag",     (DL_FUNC) &addvalue2diag,      4},
    {NULL, NULL, 0}
};

void R_init_SFSI(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
