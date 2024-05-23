#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

mat inv_qr_(mat& X);

// arma::uvec setdiff_(arma::uvec x, arma::uvec y);
// arma::vec rep_nout(vec x, int nout);
// arma::mat vec2mat(arma::vec x, int nrow, int ncol);

// arma::vec remove_elem (arma::vec X, arma::uvec ent_, int shift=1);
// arma::vec extract_elem(arma::vec X, arma::uvec ent_, int shift=1);
// arma::vec replace_elem(arma::vec X, arma::uvec ent_, arma::vec value_, int shift=1);

// arma::mat remove_rows (arma::mat X, arma::uvec row_ent_, int shift=1);
// arma::mat extract_rows(arma::mat X, arma::uvec row_ent_, int shift=1);
// arma::mat replace_rows(arma::mat X, arma::uvec row_ent_, arma::vec value, int shift=1);
  
// arma::mat replace_u_vc (arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift=1);
// arma::mat replace_uc_v (arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift=1);
// arma::mat replace_uc_vc(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift=1);
// arma::mat replace_u_v  (arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift=1);
// arma::mat replace_uv_  (arma::mat& M, arma::ivec u, arma::ivec v, arma::vec value, int shift=1);

// arma::mat extract_uv_(arma::mat& M, arma::ivec u, arma::ivec v, int shift=1);
// void as_colvec(arma::mat& M);
