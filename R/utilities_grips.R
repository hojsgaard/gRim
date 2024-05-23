#' Fast computation of covariance / correlation matrix
#'
#' @param x a numeric matrix(like object).
#' @param center,scale Should columns in x be centered and/or scaled
#' @export
fast_cov <- function(x, center=TRUE, scale=TRUE){
    A <- scale(x, center=center, scale=scale) |> crossprod()
    return(A / (nrow(x) - 1))
}

#' Genrate matrix of N(0, 1) variables
#'
#' @param n.obs Number of rows
#' @param nvar Number of columns
#' @param seed Seed for random number generator
#' 
#' @export
generate_n01 <- function(n.obs, nvar, seed=2022){
    set.seed(seed)
    d <-matrix(rnorm(n.obs * nvar), nrow=n.obs)
    d
}

solve_fun <- function(A, b=NULL){
    solve_qr <- function(A, b=NULL){
        if (is.null(b))
            solve.qr(qr(A))
        else
            solve.qr(qr(A)) %*% b
    }
    
    ## solve_ginv <- function(A, b=NULL){
    ##     if (is.null(b))
    ##         MASS::ginv(A)
    ##     else
    ##         MASS::ginv(A) %*% b
    ## }
    
    ## solveM <- function(A, b=NULL){
    ##     if (is.null(b))
    ##         Matrix::solve(A)
    ##     else
    ##         Matrix::solve(A, b)
    ## }
    
    solve_qr(A, b)
}

ggm_logL <- function(S, K, nobs){
    nobs * (log(det(K)) - sum(S * K)) / 2
}


## UTILITIES

.get.time <- function(){
    Sys.time()
}

.get.diff.time <- function(t0, units="secs"){
  unit_ <- match.arg(units, c("secs", "millisecs"))
  d <- as.numeric(difftime(Sys.time(), t0, units="secs"))
  if (identical(units, "millisecs"))
      d <- d * 1000
  d
}


matchKSigma <- function(object){
    if (!inherits(object, "gips_fit_class"))
        stop("'object' is not a gips_git_class object")
    ## sum(abs((object$K %*% object$Sigma) - diag(1, nrow(object$K))))
    a <- (object$K %*% object$Sigma) - diag(1, nrow(object$K))
    ## sqrt(sum(a^2))
    norm(a, type="F")    
}


## #' Parameter equality
## #'
## #' @param p2a,p2b Parameter lists.
## #' @param eps A small number.
## #' 
## #' @export
## #' @rdname utilities
## parmEq <- function(p2a, p2b, eps=1e-4){
##   max(abs(p2a$K - p2b$K)) / max(abs(p2a$K)) < eps && max(abs(p2a$Sigma - p2b$Sigma)) / max(abs(p2a$Sigma)) < eps
## }


#' Impose zeros in matrix entries which do not correspond to an edge.
#'
#' @param emat Edge matrix (2 x p matrix)
#' @param K Matrix; typically a concentration matrix.
#'
#' @export
impose_zero <- function(emat, K){
    emc <- as_emat_complement(emat, nrow(K))
    if (ncol(emc) == 0){
        return(K)  
    } 
    emc <- t.default(emc)
    K[rbind(emc, emc[2:1,])]  <- 0
    return(K)
}











### Internal .functions

.sol <- function(A){ ## FIXME : silly name
    .solve2x2 <- function(A){
        ##d <- A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1] ## FIXME: Handle d=0 case
        d <- A[1] * A[4] - A[2] * A[3] ## FIXME: Handle d=0 case
        #B <- A
        B <- c(0, 0, 0, 0)
        B[1] <- A[4]
        B[2] <- -A[2]
        B[3] <- -A[3]
        B[4] <- A[1]
        B <- B * (1 / d)
        dim(B) = c(2L, 2L)
        B
    }
    
    if (dim(A)[1] == 2) .solve2x2(A)
    else solve.default(A) ## FIXME: Could do chol2inv(chol(A))
}

.checkS <- function(S){
    if (!is.matrix(S)) 
        stop("Second argument is not a matrix!")
    
    if (dim(S)[1] != dim(S)[2]) 
        stop("Second argument is not a square matrix!")
}

.penalizeS <- function(S, lambda=0){
    S + diag(lambda, nrow(S))
}


## Just once
.make_rc <- function(glist, S){
    if (!is.matrix(glist))
        glist <- do.call(cbind, glist)
    vn    <- unique.default(sort(glist))
    amat  <- matrix(0L, nrow=length(vn), ncol=length(vn))
    amat[t.default(glist)] <- 1L
    amat <- amat + t.default(amat)
    dimnames(amat) <- dimnames(S)
    lt   <- lower.tri(amat)
    rrr  <- which.arr.index(lt * amat != 0)
    colnames(rrr) <- c("row", "col")
    rrr
}


#' @title Utilities for gRips
#'
#' @name utilities_grips

#' @export
#' @rdname utilities_grips
#' @param object Model object.
#' @param ... Additional arguments; currently not used.
#' @param k Penalty parameter for calculating AIC; only k=2 gives genuine AIC.
#' @param x Object to be printed.
logLik.gips_fit_class <- function(object, ...){

  nobs <- nobs(object)# unname(object$details["nobs"])
  if (is.na(nobs))
    stop("'nobs' not given; can not compute log L\n")
  
  trKS  <- unname(object$details$trKS)
  nparm <- unname(object$details$dim)
  
  ##str(list(nobs, trKS))
  out  <- c(nobs * (log(det(object$K)) - trKS) / 2)
  
  ## FIXME: Perhaps need number of parameters in specified model.
  attr(out, "nobs")  = nobs
  attr(out, "nparm") = nparm
  attr(out, "df")    = nparm
  
  class(out) <- "logLik"
  out    
}


#' @method AIC gips_fit_class
#' @export
#' @rdname utilities_grips
AIC.gips_fit_class <- function(object, ..., k=2){
    ll <- logLik(object) 
    -2 * as.numeric(ll) + k * attr(ll, "nparm")
}

#' @method BIC gips_fit_class
#' @export
#' @rdname utilities_grips
BIC.gips_fit_class <- function(object, ...){
    ll <- logLik(object) 
    -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "nparm")    
}

#' @method sigma gips_fit_class
#' @export
#' @rdname utilities_grips
sigma.gips_fit_class <- function(object, ...) {
    object$Sigma
}

#' @export
#' @rdname utilities_grips
concentration <- function(object, ...) {
    UseMethod("concentration")
}

#' @method concentration gips_fit_class
#' @export
#' @rdname utilities_grips
concentration.gips_fit_class <- function(object, ...) {
    object$K
}


#' @export
#' @rdname utilities_grips
print.gips_fit_class <- function(x, ...){
    ## cat("Method: ", x$method, "\n")
    ## xx <- c(method=x$method, eng=x$engine, x$details)
    ## print(as.data.frame(x$details))
    xx <- x$details
    xx <- as.data.frame(xx)
    xx$engine <- NULL
    xx$ncore <- NULL
    print(as.data.frame(xx))
    xx
}

#' @export
#' @rdname utilities_grips
summary.gips_fit_class <- function(object, ...) {
    xx <- as.data.frame(object$details)
    xx$engine <- NULL
    xx$ncore <- NULL
    ## xx$tpi <- with(xx, time / iter)
    ## xx$tpe <- with(xx, time / (dim-idim))
    xx
}

#' @export
#' @rdname utilities_grips 
glance.gips_fit_class <- function(x, ...) {
    as.data.frame(x$details) 
}


