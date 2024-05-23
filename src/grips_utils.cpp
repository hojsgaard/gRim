#include "RcppArmadillo.h"
// #include "general_utils.h"

using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   num_vec;
typedef Rcpp::IntegerVector   int_vec;
typedef Rcpp::CharacterVector chr_vec;

// ------------------------------------------------------------------
// General utilities - not directly related to but used by gRips.
// ------------------------------------------------------------------


//[[Rcpp::export(.c_clone)]]
SEXP clone_(SEXP& x)
{
  SEXP x2 = clone(x);
  return(x2);
}

// //[[Rcpp::export]]
// chr_vec list_names_(List lst){
//   chr_vec out(0);
//   if (lst.length() > 0)
//     out = as<chr_vec>(lst.names());
//   return out;
// }

mat as_emat2amat_(umat emat, int d){
  mat amat = zeros(d, d);
  uvec eids = sub2ind(size(amat), emat);
  vec vals = ones(emat.n_cols, 1);
  amat(eids) = vals;  
  amat = amat + amat.t();
  return amat;
}


umat as_emat_complement_(umat emat, int d){
  mat amat = as_emat2amat_(emat, d);
  amat = amat - 1;
  amat = trimatl(amat);
  amat.diag().zeros();
  uvec indices = find(amat < 0);
  umat ematc   = ind2sub(size(amat), indices);
  return ematc;  
}


// ---------------------------------------------------------------------
// *** gRips utilities ***
// ---------------------------------------------------------------------

// // [[Rcpp::export]]
// bool has_full_rank_(mat& Delta, double eps){ // FIXME
//     vec eigval;
//     mat eigvec;
    
//     eig_sym(eigval, eigvec, Delta);
//     double mev = min(eigval);
//     return (mev > eps);
//     // Rprintf("val %f sign %f mev %f\n", val, sign, mev);


// // uword rank_Delta = arma::rank(Delta); //, sqrt(datum::eps) * Delta.n_cols);
//   // return (rank_Delta >= Delta.n_cols); 	    
// }


double get_mev(mat& Delta){
    vec eigval;
    mat eigvec;
    
    eig_sym(eigval, eigvec, Delta);
    double mev = min(eigval);
    return (mev);
}


mat project_onto_G_(const mat& Delta, const umat& emc){

  mat Delta2 = Delta;
  uvec r0 = {0}, r1 = {1};
  mat emc2 = conv_to<mat>::from(emc);
  uvec s0 = conv_to<uvec>::from(emc2.rows(r0));
  uvec s1 = conv_to<uvec>::from(emc2.rows(r1));
  for (size_t j=0; j<emc.n_cols; j++){
    // Rcout << s0[j] << "  " << s1[j] << "\n"; 
    Delta2(s0[j], s1[j]) = 0;
    Delta2(s1[j], s0[j]) = 0;
  }
  return(Delta2);
}

double mnorm_one_(mat& Delta){
  rowvec s = sum(abs(Delta));
  // s.print();
  return(max(s));
}

double mnorm_maxabs_(mat& Delta){
  rowvec s = max(abs(Delta));
  // s.print();
  return(max(s));
}

mat initSigma_(mat& S)
{
  int p = S.n_rows;
  vec s = S.diag();
  mat Sigma(p, p, fill::eye);
  Sigma.diag() = s;
  return Sigma;
}

mat initK_(mat& S)
{
  int p = S.n_rows;
  vec s = S.diag();
  mat K(p, p, fill::eye);
  K.diag() = 1/s;
  return K;
}


// [[Rcpp::export]]
double ggm_logL_(mat& S, mat& K, int nobs)
{
  double trKS = accu(K % S);
  // Rf_PrintValue(wrap(trKS));
  double val, sign;
  log_det(val, sign, K);
  // double logL = nobs * (log(det(K)) - trKS) / 2;
  double logL = nobs * (val - trKS) / 2;

  return logL;
}





