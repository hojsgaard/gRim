#################################################################################
#' @title Iterative proportional fitting of graphical Gaussian model
#' 
#' @description Fit graphical Gaussian model by iterative proportional fitting.
#' 
#' @details \code{ggmfit} is based on a C implementation.  \code{ggmfitr} is
#'     implemented purely in R (and is provided mainly as a benchmark for the
#'     C-version).
#################################################################################
#' @aliases ggmfit ggmfitr
#' @param S Empirical covariance matrix
#' @param n.obs Number of observations
#' @param glist Generating class for model (a list)
#' @param start Initial value for concentration matrix
#' @param eps Convergence criterion
#' @param iter Maximum number of iterations
#' @param details Controlling the amount of output.
#' @param ... Optional arguments; currently not used
#' @return A list with \item{lrt}{Likelihood ratio statistic (-2logL)}
#'     \item{df}{Degrees of freedom} \item{logL}{log likelihood}
#'     \item{K}{Estimated concentration matrix (inverse covariance matrix)}
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{cmod}}, \code{\link{loglin}}
#' @keywords multivariate models
#' @examples
#' 
#' ## Fitting "butterfly model" to mathmark data
#' ## Notice that the output from the two fitting functions is not
#' ## entirely identical.
#' data(math)
#' glist <- list(c("al", "st", "an"), c("me", "ve", "al"))
#' d <- cov.wt(math, method="ML")
#' ggmfit (d$cov, d$n.obs, glist)
#' ggmfitr(d$cov, d$n.obs, glist)
#' 
#' @export ggmfit
ggmfit <- function(S, n.obs, glist, start=NULL, 
                   eps=1e-12, iter=1000, details=0, ...)
{
    
    glist.num <- glist
    nms_in_data <- colnames(S)

    ## Numerical (indices) representation of glist
    glist.num <- lapply(glist, match, nms_in_data)
    
    ## The used variables
    nms_in_model <- unique.default(unlist(glist))       
    
    ## Check that the used variables are in S
    zzz <- match(nms_in_model, nms_in_data)
    if (any(is.na(zzz)))
        stop("Variables ", nms_in_model[is.na(zzz)], " not in data\n")
    
    ## Get variables in the right order: The order in the S
    nms_in_model  <-  nms_in_data[sort(zzz)]
    
    ## Possibly consider only submatrix of S
    S <- S[nms_in_model, nms_in_model, drop=FALSE]
    nms_in_data <- colnames(S)

    if (is.null(start)){
        start <- diag(1/diag(S))   #print(start)
    }
        

    ## Used for calling c-code
    vn_ <- seq_along(nms_in_data)
    nvar_ <- length(vn_)

    glen_    <- sapply(glist.num, length)
    ng_      <- length(glist.num)
    
    clist.num   <- lapply(glist.num, function(x) vn_[-x])  
    clen_    <- sapply(clist.num,length)
        
    gg <- as.integer(unlist(glist.num) - 1)
    cc <- as.integer(unlist(clist.num) - 1)
        
    xxx<-.C("Cggmfit", S=S, n=as.integer(n.obs), K=start, nvar=nvar_, ngen=ng_, 
            glen=glen_, glist=gg, clen=clen_, clist=cc, 
            logL=numeric(1), eps=as.numeric(eps),
            iter=as.integer(iter), converged=as.integer(1),
            details=as.integer(details),
            PACKAGE="gRim")
    xxx <- xxx[c("logL", "K", "iter")]  
    
    dimnames(xxx$K) <- dimnames(S)
    detK  <- det(xxx$K)
    dev   <- -n.obs * log(det(S %*% xxx$K))            ## deviance to the saturated model  
    df    <-  sum(xxx$K == 0) / 2

    ## FIXME nvar_ bruges her
    out  <- list(dev=dev, df=df, detK=detK, nvar=nvar_, S=S, n.obs=n.obs)
    out   <- c(out, xxx)
    
    return(out)  
}

#' @export  
ggmfitr <- function(S, n.obs, glist, start=NULL, 
                    eps=1e-12, iter=1000, details=0, ...)
{
    ell <- function(Sigma, S, n){
        
        shdet <- function(Sigma){
            prod(eigen(Sigma)[[1]])
        }
        p <- dim(S)[1]
        const <- -n * p/2 * log(2 * pi)
        const - n/2 * log(shdet(Sigma)) - n/2 * sum(diag( solve(Sigma) %*% S )) 
    }
    
    ellK <- function (K, S, n)
    {
        value <- (n/2) * (log(det(K)) - sum(rowSums(K * S)))
        value
    }
    
    if (is.null(start)){
        K     <- diag(1/diag(S))
    } else {
        K     <- start
    }
    
    dimnames(K) <- dimnames(S)
    vn <- colnames(S); #print(vn)
    x  <- lapply(glist, match, vn)
  
    varIndex = 1:nrow(K)
    itcount = 0

    if (length(x)){
        my.complement <- function(C) return(setdiff(varIndex,C))
        x.complements <- lapply(x, my.complement)
                                        #print("x"); print(x)
                                        #print("x.comp");print(x.complements)
    
    if(length(x.complements[[1]])==0){
        return(list(K=solve(S)))
    }
        logLvec <- NULL
        repeat {    
            
      for(j in 1:length(x)){
          C     <- x[[j]]
          notC  <- x.complements[[j]]
                                        #print(C); print(S[C,C,drop=FALSE])
          K[C,C] <- solve( S[C,C,drop=FALSE] ) +
              K[C,notC,drop=FALSE]%*%solve(K[notC,notC,drop=FALSE])%*%K[notC,C,drop=FALSE]
                                        #print(K)
      }
            logL <- ellK(K,S,n.obs)
            logLvec <- c(logLvec, logL)
            itcount <- itcount+1
            if (itcount>1){
                if (logL - prevlogL < eps){  
                    converged=TRUE
                    break
                }
            } else {
                if (itcount==iter){
                    converged=FALSE
                    break
                } 
            }
            prevlogL <- logL
        }    
    }
    
    df <- sum(K[upper.tri(K)] == 0)
                                        #out <- list(K=K, logL=logL, converged=converged, itcount=itcount)
    out <- list(dev=-2*logL, df=df, logL=logL, K=K, S=S,n.obs=n.obs,
                itcount=itcount, converged=converged, logLvec=logLvec)
    return(out)
}

