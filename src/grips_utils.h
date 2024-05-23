#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   num_vec;
typedef Rcpp::IntegerVector   int_vec;
typedef Rcpp::CharacterVector chr_vec;

const int erp=0; // ERROR PRINTING
#define as_num(cc) NumericVector(cc.begin(), cc.end())

// ------------------------------------------------------------------
// General utilities - not directly related to but used by gRips.
// ------------------------------------------------------------------

SEXP clone_(SEXP& x);
// chr_vec list_names_(List lst);
mat as_emat2amat_(umat emat, int d);
umat as_emat_complement_(umat emat, int d);


// ---------------------------------------------------------------------
// *** gRips utilities ***
// ---------------------------------------------------------------------

// bool has_full_rank_(mat& Delta, double eps); // FIXME
double get_mev(mat& Delta);
  
mat project_onto_G_(const mat& Delta, const umat& emc);
double mnorm_one_(mat& Delta);
double mnorm_maxabs_(mat& Delta);
  
// int method2int_(CharacterVector method);
double ggm_logL_(mat& S, mat& K, int nobs);
mat initSigma_(mat& S);
mat initK_(mat& S);
// double get_conv_ref(const List& aux);



// -------------------------------------------------------------
// Macros common to COVIPS and CONIPS
// -------------------------------------------------------------

#define RETURN_VALUE					\
  return List::create(						\
    _["K"]     = K,						\
    _["Sigma"] = Sigma,						\
    _["logL"]  = logL,						\
    _["iter"]  = itcount,					\
    _["gap"]   = gap,						\
    _["version"] = version,					\
    _["conv_check"] = conv_check);				\
  

#define INIT_CONVERGENCE_CHECK			\
  switch (convcrit){                            \
  case 1:					\
    break;					\
  case 2:					\
    logLp = ggm_logL_(S, K, nobs);		\
    break;					\
  case 4:                                       \
    break;                                      \
  }                                             \

#define PRINT_CONV_CHECK1			                \
  if ((print > 0) && ((itcount % 1) == 0)) {			\
  logL  = ggm_logL_(S, K, nobs);				\
  sprintf(buffer, "itcount=%4i, logL=%f, nparm=%4f, mad=%f\n",	\
	  itcount, logL, nparm, mad);				\
  Rcout << buffer;						\
  }								\


#define PRINT_CONV_CHECK3						\
  if ((print > 0) && ((itcount % 1) == 0)) {				\
    sprintf(buffer, "itcount=%4i, logL=%f, logLp=%f, nparm=%4f, conv_ref=%f, conv_check=%f\n", \
	    itcount, logL, logLp, nparm, conv_ref, conv_check);		\
    Rcout << buffer;							\
  }									\

#define CONV_CHECK_LOGL_DIFF			\
  logL  = ggm_logL_(S, K, nobs);		\
  conv_check = fabs(logL - logLp)  / nparm;	\
  logLp = logL;					\

// Rprintf("logL: %14.10f logLp: %14.10f diff: %14.10f conv_check: %14.10f\n", logL, logLp, logL-logLp, conv_check); 


#define CONV_CHECK_LOGL_DIFF_REF			\
  logL  = ggm_logL_(S, K, nobs);			\
  conv_check = fabs(logL - conv_ref) / nparm;		\
  logLp = logL;						\

