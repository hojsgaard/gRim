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
extern SEXP _gRim_C_ghk2pms(void *);
extern SEXP _gRim_C_pms2ghk(void *);
extern SEXP _gRim_clone_(void *);
extern SEXP _gRim_conips_ggm_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _gRim_covips_ggm_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _gRim_fit2way_(void *, void *, void *, void *);
extern SEXP _gRim_ggm_logL_(void *, void *, void *);
extern SEXP _gRim_inv_qr_(void *);
extern SEXP _gRim_ncd_ggm_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _gRim_parm_ghk2pms_(void *);
extern SEXP _gRim_parm_normalize_ghk_(void *);
extern SEXP _gRim_RcppExport_registerCCallable(void);
extern SEXP C_ghk2pms(void *);
extern SEXP C_pms2ghk(void *);

static const R_CMethodDef CEntries[] = {
    {"Cggmfit", (DL_FUNC) &Cggmfit, 14},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_gRim_C_ghk2pms",                    (DL_FUNC) &_gRim_C_ghk2pms,                     1},
    {"_gRim_C_pms2ghk",                    (DL_FUNC) &_gRim_C_pms2ghk,                     1},
    {"_gRim_clone_",                       (DL_FUNC) &_gRim_clone_,                        1},
    {"_gRim_conips_ggm_",                  (DL_FUNC) &_gRim_conips_ggm_,                  10},
    {"_gRim_covips_ggm_",                  (DL_FUNC) &_gRim_covips_ggm_,                  10},
    {"_gRim_fit2way_",                     (DL_FUNC) &_gRim_fit2way_,                      4},
    {"_gRim_ggm_logL_",                    (DL_FUNC) &_gRim_ggm_logL_,                     3},
    {"_gRim_inv_qr_",                      (DL_FUNC) &_gRim_inv_qr_,                       1},
    {"_gRim_ncd_ggm_",                     (DL_FUNC) &_gRim_ncd_ggm_,                     10},
    {"_gRim_parm_ghk2pms_",                (DL_FUNC) &_gRim_parm_ghk2pms_,                 1},
    {"_gRim_parm_normalize_ghk_",          (DL_FUNC) &_gRim_parm_normalize_ghk_,           1},
    {"_gRim_RcppExport_registerCCallable", (DL_FUNC) &_gRim_RcppExport_registerCCallable,  0},
    {"C_ghk2pms",                          (DL_FUNC) &C_ghk2pms,                           1},
    {"C_pms2ghk",                          (DL_FUNC) &C_pms2ghk,                           1},
    {NULL, NULL, 0}
};

void R_init_gRim(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
