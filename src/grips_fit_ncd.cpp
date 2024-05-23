#include "RcppArmadillo.h"
#include "grips_utils.h"
// #include "general_utils.h"
#include "grips_arma_utils.h"
// #include "convergence.h"
#include "grips_precision.h"
#include <iostream>
#include <vector>
#include <algorithm> // for std::nth_element
#include <stdio.h>

using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   num_vec;
typedef Rcpp::IntegerVector   int_vec;
typedef Rcpp::CharacterVector chr_vec;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))

// ### ###################################################
// ### UTILITY FUNCTIONS FOR NCD
// ### ###################################################

double duality_gap_(mat& S, mat& K, int nobs){
  mat KS = K * S; 
 
  double val, sign;  
 log_det(val, sign, KS); 

  // double out = nobs * (accu(K % S) - log(det(KS)) - S.n_rows) / 2; // FIXME Fragile  
     double out = nobs * ( - val + accu(K % S) - S.n_rows) / 2; //  Fragile 
	
     out = fabs(out); // to avoid negative numbers by rounding error
 
  return out;
}


// ### ###################################################
// ### NCD FUNCTIONS
// ### ###################################################

// ### ###################################################
// ### update rows of Sigma and K
// ### ###################################################

// *** Used by outerloop1 and outerloop2

void update_Sigma_row_(int u, mat& Sigma, const mat& amat, int nobs, int print=0){

  // FIXME Need not be computed each time...
  uvec u_    = {(unsigned int) u};      // convert int to uvec
  uvec ub_   = find(amat.rows(u_) > 0); // Returns column vector
  int  deg_u = accu(amat.rows(u_));

  if (print >= 4){
    Rprintf(">>>> Updating Sigma for u=%i with degree %i\n", u, deg_u);
  }

  vec w(Sigma.n_cols, fill::zeros);  
  vec beta_star;
  
  if (ub_.n_rows > 0){  
    mat AA    = Sigma.submat(ub_, ub_);
    mat s_u   = Sigma.cols(u_);
    vec s_ubu = s_u.rows(ub_);

    // if (deg_u > nobs - 1){
    //   beta_star = pinv(AA) * s_ubu;       // Rprintf("using pinv\n");
    // } else {
    //   beta_star = solve(AA, s_ubu);       // Rprintf("using solve\n");      
    // }

    beta_star = solve(AA, s_ubu);       // Rprintf("using solve\n");
      
    w = Sigma.cols(ub_) * beta_star;
  }

  // Rprintf("inserting:\n");
  double sigma_uu = Sigma(u, u); // Store this because element is overwritten below    
  Sigma.col(u) = w;
  Sigma.row(u) = w.t();
  Sigma(u, u)  = sigma_uu;       // Restore element
}


// *** Used by outerloop2

bool shall_update(int u, mat& K, mat& amat, double eps=0.01){
  uvec u_    = {(unsigned int) u};    
  uvec ur_   = find(amat.rows(u_) == 0); // Returns column vector  
  uvec locate_u = find(ur_ == u);
  ur_.shed_rows(locate_u);

  mat K_uru      = K.submat(ur_, u_);    // Column vector
  double mno = mnorm_one_(K_uru);  
  return (mno > eps);
}

void update_K_row_(int u, mat& Sigma, mat& K, const mat& amat, int print=0){

    uvec u_    = {(unsigned int) u};    
    uvec ub_   = find(amat.rows(u_) > 0); // Returns column vector      
    uvec uc_   = arma::linspace<arma::uvec>(0, K.n_cols - 1, K.n_cols);
    uc_.shed_rows(u_);
    int  deg_u = accu(amat.rows(u_));
    if (print >= 4){
      Rprintf(">>>> Updating K for u=%i with degree %i\n", u, deg_u);
      Rprintf(">>>> ub_: ");  ub_.t().print();
    }
    
    double k_uu     = as_scalar(K(u, u));    
    double sigma_uu = as_scalar(Sigma(u, u));
    
    mat K_ucu      = K.submat(uc_, u_);
    mat K_ucuc     = K.submat(uc_, uc_);
    mat Sigma_ucu  = Sigma.submat(uc_, u_);
    mat CC2        = K_ucuc - K_ucu * (trans(K_ucu) / k_uu);
    mat DD2        = CC2 * Sigma_ucu;
    mat k2_uu      = 1 / (sigma_uu - trans(Sigma_ucu) * DD2);
    mat K2_ucu     = as_scalar(k2_uu) * DD2;
    
    // DO UPDATE    
    mat RR2, RR, K2_ucuc;
    
    K(u, u)           = as_scalar(k2_uu);
    K.submat(uc_, u_) = -K2_ucu;
    K.submat(u_, uc_) = trans(-K2_ucu);
    
    RR2 = K2_ucu * (trans(K2_ucu) / as_scalar(k2_uu));
    RR  = K_ucu  * (trans(K_ucu)  / as_scalar(k_uu));
    
    K2_ucuc  = K_ucuc + RR2 - RR;
    K.submat(uc_, uc_) = K2_ucuc;    
}



// ### ###################################################
// ### OUTERLOOP1
// ### ###################################################


// Update all nodes once
void ncd_inner1_update_Sigma_(mat& Sigma, mat& amat, int nobs, int print=0){

  if (print >= 4){
    Rprintf(">>>> Running ncd_inner1_update_Sigma\n");
  }
  
  for (size_t u=0; u<amat.n_rows; u++){
    update_Sigma_row_(u=u, Sigma=Sigma, amat=amat, nobs=nobs, print=print);
  }
}

List ncd_outer1_(mat& Sigma, mat& K, umat& emat, umat& emat_c, mat& amat,
		 int& nobs, double& eps, int max_visits, int& n_visits, int print=0){

  if (print >=2){
    Rprintf(">> Running ncd_outer1\n");
  }
  
  double conv_check, mno;
  mat Sigma_prev = diagmat(Sigma.diag());
  int n_vars = Sigma.n_rows, count=0, n_upd = n_vars;
  
  while (true) {
    ncd_inner1_update_Sigma_(Sigma=Sigma, amat=amat, nobs=nobs, print=print);
    mat Delta = Sigma - Sigma_prev;
    mno       = mnorm_one_(Delta);
    n_visits += n_upd;
    count++;

    if (print >=3){
      Rprintf(">>> ncd_outer1 count: %4d max_visits: %7d n_visits: %7d n_upd: %5d mno: %10.6f eps: %8.6f\n",
	      count, max_visits, n_visits, n_upd, mno, eps);
    }
    
    Sigma_prev = Sigma;
    conv_check = mno;
    if ((n_visits == max_visits) || (conv_check < eps)) break;
  }
  return List::create(_["iter"]  = n_visits, _["mad"]=mno); //FIXME mad should be conv_crit		
}


// ### ###################################################
// ### OUTERLOOP2
// ### ###################################################

void ncd_inner2_update_Sigma_K_(mat& Sigma, mat& K, mat& amat, int nobs,
				int &n_upd, double eps=0.01, int print=0){
  if (print >= 4){
    Rprintf(">>>> Running ncd_inner2_update_Sigma_K\n");
   }

  for (size_t u=0; u<amat.n_rows; u++){
    if (shall_update(u=u, K=K, amat=amat, eps=eps)){
      n_upd++;
      update_Sigma_row_(u=u, Sigma=Sigma,      amat=amat, nobs=nobs, print=print);    
      update_K_row_    (u=u, Sigma=Sigma, K=K, amat=amat, print=print);
    }   
  }
}

List ncd_outer2_(mat& Sigma, mat& K, umat& emat, umat& emat_c, mat& amat, int nobs, double& eps, 
		 int max_visits, int& n_visits, 
		 int& n_upd, int print=0){

  if (print >=2){
    Rprintf(">> Running ncd_outer2\n");
  }

  int count=0;
  double mno, conv_crit;
  while (true){
    n_upd = 0;
    ncd_inner2_update_Sigma_K_(Sigma=Sigma, K=K, amat=amat, nobs=nobs,
			       n_upd=n_upd, eps=eps, print=print);

    n_visits += n_upd;
    mat Delta = K - project_onto_G_(K, emat_c);
    mno = mnorm_one_(Delta);
    conv_crit = mno;
    count++;
    
    if (print >=3){
      // double val, sign;
      // log_det(val, sign, Sigma);
      
      Rprintf(">>> ncd_outer2 count: %4d max_visits: %7d n_visits: %7d n_upd: %5d mno: %10.6f eps: %8.6f\n",
	      count, max_visits, n_visits, n_upd, mno, eps);
    }
    
    if ((n_visits == max_visits) || (conv_crit < eps)) break;    
  }
  return List::create(_["iter"] = n_visits, _["conv_crit"] = conv_crit);		
}


// ### ###################################################
// ### MAIN NCD FUNCTION 
// ### ###################################################

//[[Rcpp::export(.c_ncd_ggm_)]]
List ncd_ggm_(mat& S, List& elst, umat& emat, int& nobs,
	      mat K,       
	      int maxit, double& eps, int& convcrit, int print, List& aux){
  
  int version      = aux["version"];
  bool converged;
  int n_vars = S.n_cols;
  int max_visits = n_vars * maxit, n_visits = 0;
  
  umat emat_c = as_emat_complement_(emat-1, n_vars);
  // mat amat    = as_emat2amat_(emat-1, n_vars);
  mat amat = aux["amat"];
  // amat.print();
  mat Sigma   = S, K2, Delta;
  List res1, res2;
  double logL, gap=-1.0, conv_check, eps2, mno;
  int iter1, iter2, itcount, n_upd=0;

  double eps1 = 2 * eps / nobs;
  eps2 = MIN(eps1, 1.0/Sigma.n_rows);  
  
  switch (version){
  case 0:
    res1 = ncd_outer1_(Sigma=Sigma, K=K, emat=emat, emat_c=emat_c, amat=amat,
		       nobs=nobs, eps=eps1, max_visits=max_visits, n_visits=n_visits, print=print);
    break;
  case 1:
    res1 = ncd_outer1_(Sigma=Sigma, K=K, emat=emat, emat_c=emat_c, amat=amat,
		       nobs=nobs, eps=eps2, max_visits=max_visits, n_visits=n_visits, print=print);
    break;
  }
  
  iter1 = res1["iter"];
  
  if (print>=2)
    Rprintf(">> ncd_outer1 visits : %d\n", iter1);

  // FIXME mev
  double mev = get_mev(Sigma);
  // Rprintf("mev: %f\n", mev);
  if (mev > eps2){
    converged = true;
    K = inv_qr_(Sigma);
  } else {
    converged = false;
    REprintf("NCD not converged: Smallest eigenvalue = %9.6f\n", mev);	
    itcount = iter1;
  }									
  
  if (converged){
    switch (version){      
    case 0:     
      K     = inv_qr_(Sigma);
      K2    = project_onto_G_(K, emat_c);
      Delta = K - K2;
      mno   = mnorm_one_(Delta);
      
      if (print>=3)
	Rprintf(">>> sncd mno : %14.10f\n", mno);

      logL       = ggm_logL_(S, K, nobs);
      conv_check = mno;      
      gap        = -1; 
      // K    = K2;  // NOTE K2 is not returned...      
      itcount = iter1 + 1;
      break;
      
    case 1: // FULL VERSION
      K = inv_qr_(Sigma);
      n_visits = 0;
      res2 = ncd_outer2_(Sigma=Sigma, K=K, emat=emat, emat_c=emat_c, amat=amat, nobs=nobs, eps=eps2,
			 max_visits=max_visits, n_visits=n_visits, 
			 n_upd=n_upd, print=print);
      iter2 = res2["iter"];
      if (print>=2)
	Rprintf(">> ncd_outer2 visits: %d\n", iter2);
      
      K2    = project_onto_G_(K, emat_c);
      Delta = K - K2;
      mno   = mnorm_one_(Delta);
      if (print>=3)
	Rprintf(">>> ncd mno : %14.10f\n", mno);
      conv_check = mno;

      if (iter2 <= max_visits){ // Then K is posdef	
	logL = ggm_logL_(S, K2, nobs);
	gap  = duality_gap_(Sigma, K2, nobs);
	K    = K2;
      } else {
	REprintf("Algorithm may not have converged\n");
	// K = NA; upper_limit_logL = formel (23)// FIXME
      }
      itcount = iter1 + iter2;  
      break;
      
    default:
      Rprintf("'version' must be 0, 1\n");  
    }
  } else {
    logL       = -1;
    conv_check = -1;
  }


  itcount = itcount / n_vars;
  
  return List::create(						
    _["K"]     = K,						
    _["Sigma"] = Sigma,						
    _["logL"]  = logL,						
    _["iter"]  = itcount,					
    _["gap"]   = gap,
    _["mev"]   = mev,
    _["version"]    = version,
    _["converged"]  = converged, 
    _["conv_check"] = conv_check);				
  
}

