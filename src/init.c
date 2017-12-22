#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void Cggmfit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP C_ghk2pms(SEXP);
extern SEXP C_pms2ghk(SEXP);
extern SEXP gRim_ghk2pmsParms_(SEXP);
extern SEXP gRim_normalize_ghkParms_(SEXP);
extern SEXP gRim_pms2ghkParms_(SEXP);
extern SEXP gRim_updateA(SEXP, SEXP, SEXP, SEXP);
extern SEXP gRim_update_ghkParms_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"Cggmfit", (DL_FUNC) &Cggmfit, 14},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"C_ghk2pms",                (DL_FUNC) &C_ghk2pms,                1},
    {"C_pms2ghk",                (DL_FUNC) &C_pms2ghk,                1},
    {"gRim_ghk2pmsParms_",       (DL_FUNC) &gRim_ghk2pmsParms_,       1},
    {"gRim_normalize_ghkParms_", (DL_FUNC) &gRim_normalize_ghkParms_, 1},
    {"gRim_pms2ghkParms_",       (DL_FUNC) &gRim_pms2ghkParms_,       1},
    {"gRim_updateA",             (DL_FUNC) &gRim_updateA,             4},
    {"gRim_update_ghkParms_",    (DL_FUNC) &gRim_update_ghkParms_,    9},
    {NULL, NULL, 0}
};

void R_init_gRim(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
