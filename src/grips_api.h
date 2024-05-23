#include "RcppArmadillo.h"
#include <iostream>
#include <vector>
#include <algorithm> // for std::nth_element

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::interfaces(r,cpp)]]

using namespace Rcpp;
using namespace arma;


// arma::mat sol2x2_(arma::mat& A){
//   arma::mat B(2, 2); 
//   double d  = A(0,0) * A(1,1) - A(0,1) * A(1, 0);
//   double di = 1/d;
//   B(0, 0) =  A(1, 1) * di;
//   B(1, 1) =  A(0, 0) * di;
//   B(0, 1) = -A(0, 1) * di;
//   B(1, 0) = -A(1, 0) * di;
//   return(B);	
// }


// //[[Rcpp::export(.c_sol)]]
// arma::mat sol_(arma::mat& A){
//   if (A.n_rows == 2) return sol2x2_(A);
//   else return inv(A);
// }


// //[[Rcpp::export(.c_getdlogL)]]
// double getdlogL_(const arma::mat Scc, const arma::mat& Sigcc, 
// 		  const arma::mat& Kt, const arma::mat& Ks){
//   double v1 = log(det(Kt * Sigcc));
//   //double v1 = log_det(Kt * Sigcc); // FIXME should work...
//   double v2 = 0;
//   for (int i=0; i<=1; ++i)
//     for (int j=0; j <= 1; ++j)
//       v2 += (Ks(i,j) - Kt(i,j)) * Scc(i,j);
//   return v1 + v2;
// }


// ******************************************
//
// get-functions
//
// ******************************************

double get_suv_adjust_(arma::mat Scc, double luv){
  return (2 * Scc(0, 0) * Scc(1, 1) * luv) / (sqrt(1 + 4 * luv * luv * Scc(0, 0) * Scc(1, 1)) + 1); 
}

inline arma::mat getKstar_(arma::mat& Sigmacc){
  return inv(Sigmacc); //   return sol_(Sigmacc);
}

inline arma::mat getKtilde_(arma::mat& Saux){
  return inv(Saux); //return sol_(Saux);
}

inline arma::mat getLaux_(arma::mat& Kcc, arma::mat& Kstar){
  arma::mat out = Kcc - Kstar;	
  return out;
}

inline arma::mat getKupd_(arma::mat& Ktilde, arma::mat& Laux){
  arma::mat out = Ktilde + Laux;	
  return out;
}

inline arma::mat getHaux_(arma::mat& Saux, arma::mat& Kstar){
  arma::mat out = Kstar - Kstar * Saux * Kstar;
  return out;
}

inline arma::mat getSaux_ips_(arma::mat& Scc, arma::mat& Laux){
  //Rcout << "getSaux_ips_" << std::endl;
  return Scc;
}

inline arma::mat getSaux_lasso_(arma::mat& Scc, arma::mat& Laux, arma::mat& lambdacc){
  //Rcout << "getSaux_lasso_" << std::endl;
  double suv = Scc(0, 1);
  double luv = Laux(0, 1);
  double lambdauv = lambdacc(0, 1);
  double suv_star = get_suv_adjust_(Scc + lambdacc, luv);

  double suv_tilde = suv_star;
  if (suv + lambdauv < suv_star){
    suv_tilde = suv + lambdauv;
  } else if (suv - lambdauv > suv_star){
    suv_tilde = suv - lambdauv;
  }

  arma::mat Saux = Scc + lambdacc; // FIXME: Can be avoided by
				   // adjusting diag of S once and for
				   // all
  Saux(1, 0) = suv_tilde;
  Saux(0, 1) = suv_tilde;  
  return Saux; 
}

inline arma::mat getSaux_mtp2_(arma::mat& Scc, arma::mat& Laux){
  //Rcout << "getSaux_mtp2_" << std::endl;
  double suv = Scc(0, 1);
  double luv = Laux(0, 1);

  double limuv = suv / (Scc(0, 0) * Scc(1, 1) - suv * suv);
  double suv_tilde = 0;

  if (luv <= limuv){
    suv_tilde = suv;
    //Rcout << "case 1 : " << suv_tilde << std::endl;
  } else {
    suv_tilde = get_suv_adjust_(Scc, luv);
    //Rcout << "case 2 : " << suv_tilde << std::endl;	
  }

  arma::mat Saux = Scc;
  //Rcout << "before : " << Saux(1, 0) << std::endl;
  Saux(1, 0) = suv_tilde;
  Saux(0, 1) = suv_tilde;
  //Rcout << Saux << std::endl;
  return Saux; 
}

inline arma::mat getSaux_hybrid_(arma::mat& Scc, arma::mat& Laux, arma::mat& lambdacc, bool pos=false){
  //Rcout << "getSaux_mtp2_" << std::endl;
  double suv = Scc(0, 1);
  double luv = Laux(0, 1);
  double lambdauv = lambdacc(0, 1);
  
  double suv_star = get_suv_adjust_(Scc + lambdacc, luv);
  double suv_tilde = suv_star;

  if (suv + lambdauv < suv_star){
    suv_tilde = suv + lambdauv;
  } else if (suv - lambdauv > suv_star){
    suv_tilde = suv - lambdauv;
  }

  arma::mat Saux = Scc + lambdacc;

  if (pos){
    double limuv = suv_tilde / (Scc(0, 0) * Scc(1, 1) - suv * suv);
    if (luv <= limuv){
      ; //suv_tilde = suv_tilde;
    } else {
      suv_tilde = suv_star;
    }
  }
  Saux(1, 0) = suv_tilde;
  Saux(0, 1) = suv_tilde;
  //Rcout << Saux << std::endl;
  return Saux; 
}


inline arma::mat getSaux_(arma::mat& Scc, arma::mat& Laux, arma::mat& lambdacc,
		   int imeth=0,
		   bool pos=false){
  
  switch(imeth){
  case 0 :    //Rcout << "ips" << std::endl;
    return getSaux_ips_(Scc, Laux); break;	
  case 1:     //Rcout << "mtp2" << std::endl;    
    return getSaux_mtp2_(Scc, Laux); break;
  case 2:     //Rcout << "lasso" << std::endl;    
    return getSaux_lasso_(Scc, Laux, lambdacc); break;
  case 3:     //Rcout << "hybrid" << std::endl;    
    return getSaux_hybrid_(Scc, Laux, lambdacc, pos); break;    
  default:
    Rcout << "getSaux - dunnowhattodo" << std::endl;    
  }
  //Rcout << "what...";
  return Scc; // OK; Never reached 
}
