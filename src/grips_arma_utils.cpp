#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
mat inv_qr_(mat& X){
  arma::mat Q, R;
  arma::qr(Q,R,X);
  mat Ri = inv(R);
  mat out = Ri * Q.t();
  return(out);
}

// // [[Rcpp::export]]
// arma::uvec setdiff_(arma::uvec x, arma::uvec y){    
//     // x = arma::unique(x);
//     // y = arma::unique(y);    
//     for (size_t j = 0; j < y.n_elem; j++) {
//         arma::uvec q1 = arma::find(x == y[j]);
//         if (!q1.empty()) {
//             x.shed_row(q1(0));
//         }
//     }
//     return x;
// }

// // [[Rcpp::export]]
// arma::vec rep_nout(vec x, unsigned int nout){
//   if (x.n_elem == nout){
//     return(x);
//   }
//   vec x_(nout);
//   size_t k = 0;
//   for (size_t i=0; i<nout; i++){
//     x_(i) = x(k++);
//     if (k == x.n_elem)
//       k = 0;
//   }
//   return(x_);    
// }

// // [[Rcpp::export]]
// arma::mat vec2mat(arma::vec x, int nrow, int ncol) {
//     arma::vec x_ = rep_nout(x, nrow * ncol);
//     arma::mat y(x_);
//     y.set_size(nrow, ncol);
//     return y;
// }

// // [[Rcpp::export]]
// arma::vec remove_elem(arma::vec X, arma::uvec ent_, int shift=1) {
//   X.shed_rows(ent_ - shift); // remove elements
//   return(X);
// }

// // [[Rcpp::export]]
// arma::vec extract_elem(arma::vec X, arma::uvec ent_, int shift=1) {
//     return X.rows(ent_ - shift);
// }

// // [[Rcpp::export]]
// arma::vec replace_elem(arma::vec X, arma::uvec ent_, arma::vec value_, int shift=1) {
//     X.rows(ent_ - shift) = value_;
//     return(X);
// }

// // [[Rcpp::export]]
// arma::mat remove_rows(arma::mat X, arma::uvec row_ent_, int shift=1) {
//     X.shed_rows(row_ent_ - shift); // remove rows
//     return(X);
// }

// // [[Rcpp::export]]
// arma::mat extract_rows(arma::mat X, arma::uvec row_ent_, int shift=1) {
//     return (X.rows(row_ent_ - shift));
// }

// // [[Rcpp::export]]
// arma::mat replace_rows(arma::mat X, arma::uvec row_ent_, arma::vec value, int shift=1){
//     arma::mat value_ = vec2mat(value, row_ent_.n_elem, X.n_cols);    
//     X.rows(row_ent_ - shift) = value_;
//     return(X);
// }

// // [[Rcpp::export]]
// arma::mat replace_u_vc(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift=1){
//     arma::uvec vc_  = arma::linspace<arma::uvec>(0, M.n_cols - 1, M.n_cols);
//     arma::uvec vc  = setdiff_(vc_, v - shift);
//     arma::mat value_ = vec2mat(value, u.n_elem, vc.n_elem);    
//     // Rprintf("u, uc:\n"); u.t().print(), uc.t().print();
//     M.submat(u - shift, vc) = value_;
//     return(M);
// }

// // M[-u, u] <- value
// // [[Rcpp::export]]
// arma::mat replace_uc_v(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift=1){

//   arma::uvec vr_  = arma::linspace<arma::uvec>(0, M.n_rows - 1, M.n_rows);
//   arma::uvec uc  = setdiff_(vr_, u - shift); 
//   arma::mat value_ = vec2mat(value, uc.n_elem, v.n_elem);
//   // Rprintf("u, uc:\n"); u.t().print(); uc.t().print(); v.t().print();
//   // M.print();
//   // value_.print();

//   M.submat(uc, v - shift) = value_; 
//   return(M);
// }

// // [[Rcpp::export]]
// arma::mat replace_uc_vc(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift=1){
//   arma::uvec vr_  = arma::linspace<arma::uvec>(0, M.n_rows - 1, M.n_rows);
//   arma::uvec vc_  = arma::linspace<arma::uvec>(0, M.n_cols - 1, M.n_cols);
//   arma::uvec uc  = setdiff_(vr_, u - shift);
//   arma::uvec vc  = setdiff_(vc_, v - shift);
//   arma::mat value_ = vec2mat(value, uc.n_elem, vc.n_elem);
//   // Rprintf("u, uc:\n"); u.t().print(), uc.t().print();
//   M.submat(uc, vc) = value_; 
//   return(M);
// }

// // [[Rcpp::export]]
// arma::mat replace_u_v(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift=1){
//     arma::mat value_ = vec2mat(value, u.n_elem, v.n_elem);    
//     // Rprintf("u, uc:\n"); u.t().print(), uc.t().print();
//     M.submat(u - shift, v - shift) = value_;
//     return(M);
// }



// // [[Rcpp::export]]
// arma::mat replace_uv_(arma::mat& M, arma::ivec u, arma::ivec v, arma::vec value, int shift=1){

//   uvec uu, vv;
//   int u_code = 1, v_code=1;
//   if (u(0) < 0){
//     u_code = 0;
//   } 

//   if (v(0) < 0){
//     v_code = 0;
//   } 

//   uu = conv_to<uvec>::from(abs(u));
//   vv = conv_to<uvec>::from(abs(v));      
//   // Rprintf("uu:\n"); u.print(); uu.print();
//   // Rprintf("vv:\n"); v.print(); vv.print();
//   // Rcout << "u_code: " << u_code << " v_code: " << v_code << "\n";
//   // value.print();
    
//   // Rcout << "u_code: " << u_code << " v_code: " << v_code << "\n";  
//   if ((u_code == 1) && (v_code == 1))
//     return(replace_u_v(M, uu, vv, value, shift=shift));

//   if ((u_code == 0) && (v_code == 1))
//     return(replace_uc_v(M, uu, vv, value, shift=shift));

//   if ((u_code == 1) && (v_code == 0))
//     return(replace_u_vc(M, uu, vv, value, shift=shift));

//   if ((u_code == 0) && (v_code == 0))
//     return(replace_uc_vc(M, uu, vv, value, shift=shift));
//   // Never reached
//   return(M);
// }

// // [[Rcpp::export]]
// void as_colvec(arma::mat& M){
//   M.set_size(M.n_rows * M.n_cols, 1);
// }

// // [[Rcpp::export]]
// arma::mat extract_uv_(arma::mat& M, arma::ivec u, arma::ivec v, int shift=1){

//   uvec uu, vv;
//   int u_code = 1, v_code=1;

//   if (u(0) < 0){
//     u_code = 0;
//   }
//   if (v(0) < 0){
//     v_code = 0;
//   }
  
//   uu = conv_to<uvec>::from(abs(u));
//   vv = conv_to<uvec>::from(abs(v));
      
//   int NA_num = -2147483648;

//   // Rprintf("uu:\n"); u.print(); uu.print(); 
//   // Rprintf("vv:\n"); v.print(); vv.print(); 
//   // Rcout << "u_code: " << u_code << " v_code: " << v_code << "\n";
  
//   if ((u_code == 1) && (v_code == 1)){
//     return(M.submat(uu - shift, vv - shift));
//   }

//   mat out = M;  
//   if ((u_code == 0) && (v_code == 1)){

//     out = out.cols(vv - shift);
    
//     if (u(0) != NA_num){     
//       out.shed_rows(uu - shift);
//     }
//     return(out);
//   }
  
//   if ((u_code == 1) && (v_code == 0)){
//     out = out.rows(uu - shift);
//     if (v(0) != NA_num){     
//       out.shed_cols(vv - shift);
//     }
//     return(out);
//   }

//   if ((u_code == 0) && (v_code == 0)){
//     if (v(0) != NA_num){     
//       out.shed_cols(vv - shift);
//     }
//     if (u(0) != NA_num){     
//       out.shed_rows(uu - shift);
//     }
//     return(out);
//   }

//   // Never reached
//   return(M);
// }

