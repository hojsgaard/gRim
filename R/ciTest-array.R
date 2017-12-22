## ########################################################
##
## CIP test in a table
## <x>  : table
##
## ########################################################

## FIXME: ciTest-array:
## FIXME: Clean up code;
## FIXME: use modern tabXXX functions;
## FIXME: add parametric bootstrap

#' @title Test for conditional independence in a contingency table
#' 
#' @description Test for conditional independence in a contingency table
#'     represented as an array.
#'
#' @name ciTest-array
#' 
#' @details \code{set} can be 1) a vector or 2) a right-hand sided formula in
#'     which variables are separated by '+'. In either case, it is tested if the
#'     first two variables in the \code{set} are conditionally independent given
#'     the remaining variables in \code{set}.  (Notice an abuse of the '+'
#'     operator in the right-hand sided formula: The order of the variables does
#'     matter.)
#' 
#' If \code{set} is \code{NULL} then it is tested whether the first two
#' variables are conditionally independent given the remaining variables.
#' 
#' @param x An array of counts with named dimnames.
#' @param set A specification of the test to be made. The tests are of the form u
#'     and v are independent condionally on S where u and v are variables and S
#'     is a set of variables. See 'details' for details about specification of
#'     \code{set}.
#' @param statistic Possible choices of the test statistic are \code{"dev"} for
#'     deviance and \code{"X2"} for Pearsons X2 statistic.
#' @param method Method of evaluating the test statistic. Possible choices are
#'     \code{"chisq"}, \code{"mc"} (for Monte Carlo) and \code{"smc"} for
#'     sequential Monte Carlo.
#' 
#' @param adjust.df Logical. Should degrees of freedom be adjusted for sparsity?
#' @param slice.info Logical. Should slice info be stored in the output?
#' @param L Number of extreme cases as stop criterion if method is \code{"smc"}
#'     (sequential Monte Carlo test); ignored otherwise.
#' @param B Number (maximum) of simulations to make if method is \code{"mc"} or
#'     \code{"smc"} (Monte Carlo test or sequential Monte Carlo test); ignored
#'     otherwise.
#' @param ...  Additional arguments.
#' @return An object of class 'citest' (which is a list).
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{ciTest}}, \code{\link{ciTest.data.frame}},
#'     \code{\link{ciTest_df}}, \code{\link{ciTest.list}},
#'     \code{\link{ciTest_mvn}}, \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#'
#' data(lizard)
#'
#' ## lizard is has named dimnames
#' names( dimnames( lizard ))
#' ## checked with
#' is.named.array( lizard )
#' 
#' ## Testing for conditional independence:
#' # the following are all equivalent:
#' ciTest(lizard, set=~diam + height + species)
#' # ciTest(lizard, set=c("diam", "height", "species"))
#' # ciTest(lizard, set=1:3)
#' # ciTest(lizard)
#' # (The latter because the names in lizard are as given above.)
#'
#' ## Testing for marginal independence
#' ciTest(lizard, set=~diam + height)
#' ciTest(lizard, set=1:2)
#'
#' ## Getting slice information:
#' ciTest(lizard, set=c("diam", "height", "species"), slice.info=TRUE)$slice
#' 
#' ## Do Monte Carlo test instead of usual likelihood ratio test. Different
#' # options:
#'
#' # 1) Do B*10 simulations divided equally over each slice: 
#' ciTest(lizard, set=c("diam", "height", "species"), method="mc", B=400)
#' # 2) Do at most B*10 simulations divided equally over each slice, but stop
#' # when at most L extreme values are found
#' ciTest(lizard, set=c("diam", "height", "species"), method="smc", B=400)
#' 
#' @rdname ciTest-array
ciTest_table <- function(x, set=NULL, statistic="dev", method="chisq", adjust.df=TRUE, slice.info=TRUE, L=20, B=200, ...){

  statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))
  method    <- match.arg(toupper(method),    c("CHISQ", "MC", "SMC"))

  if (is.null(set)){
    set <- names( dimnames(x) )
  } else {
    if ( inherits(set, "integer") || inherits(set, "numeric") ){
      x   <- tableMargin(x, set)
    } else
    if (inherits(set,c("formula","character"))){
      set <- unlist(rhsFormula2list(set))
      vn  <- names(dimnames(x))
      set <- vn[pmatch(set, vn)]
      x   <- tableMargin(x, set)
    }
  }

  switch(method,
         "CHISQ"={
           .CI_X2_prim(x, statistic=statistic, adjust.df=adjust.df, slice.info=slice.info)
         },
         "MC"=,"SMC"={
           .CI_SMC_prim(x, statistic=statistic, method=method, slice.info=slice.info, L=L, B=B)
         }
         )
}

###
### CIP test; asymptotic, based on either deviance or Pearsons X2
###

.CI_X2_prim <- function(x, statistic="DEV", adjust.df=TRUE, slice.info=TRUE){

  statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))

  vn     <- names(dimnames(x))
  dn     <- dim(x)
  v1     <- vn[1]
  v2     <- vn[2]
  R      <- vn[-(1:2)]
  dim.v1 <- dn[1]
  dim.v2 <- dn[2]
  dim.R  <- prod(dn[-(1:2)])

  t.v1R <- tableMargin(x, c(v1, R))
  t.v2R <- tableMargin(x, c(v2, R))

  ## Fit table
  if (length(R)){
    t.R   <- tableMargin(x, R)
    fit.table <- tablePerm(tableOp(tableOp(t.v1R, t.v2R), t.R, "/"), vn)
  } else {
    fit.table <- tablePerm(tableOp(t.v1R, t.v2R), vn)/sum(x)
  }

    ## Evaluate test statistic
    ## FIXME There are functions for that in other functions
    if (statistic=="DEV"){             ## Deviance
        tobs  <- 2 * x * log(x / fit.table)
    } else {                           ## Pearson X2
        tobs <- (x - fit.table)^2 / fit.table
    }
    tobs[!is.finite(tobs)] <- 0
    tobsGlobal <- sum(tobs)
    
    ## Calculate df with or without adjustment for sparsity
    if (adjust.df){
        t.v1Rmat    <- matrix(t.v1R, nrow=dim.R, byrow=TRUE)
        t.v2Rmat    <- matrix(t.v2R, nrow=dim.R, byrow=TRUE)

        zzz <- (t.v1Rmat > 0) * 1
        dim.v1.adj <- if (!is.null(dim(zzz))) rowSums(zzz) else sum(zzz)
                    
        zzz <- (t.v2Rmat > 0) * 1
        dim.v2.adj <- if (!is.null(dim(zzz))) rowSums(zzz) else sum(zzz)
        
        d1         <- dim.v1.adj - 1
        d1[d1 < 0] <- 0
        d2         <- dim.v2.adj - 1
        d2[d2 < 0] <- 0
        dofSlice   <- d1 * d2
    } else {
        dofSlice <- rep.int((dim.v1 - 1) * (dim.v2 - 1), dim.R)
    }
    dofGlobal <- sum(dofSlice)
    pGlobal   <- 1 - pchisq(tobsGlobal, dofGlobal)
    
    tobsSlice <- rowSumsPrim(matrix(tobs, nrow=dim.R, byrow=TRUE))
    pSlice    <- 1 - pchisq(tobsSlice, df=dofSlice)
    
    if (length(R) && slice.info){
        sliceInfo <- list(statistic=tobsSlice, p.value=pSlice, df=dofSlice)
        des   <- expand.grid(dimnames(x)[-(1:2)])
        slice <- cbind(as.data.frame(sliceInfo[1:3]), des)
    } else {
        slice <- NULL
    }
    ans <- list(statistic=tobsGlobal, p.value=pGlobal, df=dofGlobal, statname=statistic,
                method="CHISQ",
                adjust.df=adjust.df, varNames=names(dimnames(x)),slice=slice)
    class(ans) <- "citest"
    ans
}


##        if (!is.null(dim(zzz)))
##            dim.v1.adj  <- rowSums(zzz)
##        else
##            dim.v1.adj  <- sum(zzz)
        ##
        
##        if (!is.null(dim(zzz)))
##            dim.v2.adj  <- rowSums(zzz)
##        else
##            dim.v2.adj  <- sum(zzz)
##        


###
### CIP test; exact, based on sequential monte carlo
###

.CI_SMC_prim <- function(x, statistic="DEV", method="SMC", L=50, B=200, slice.info=FALSE){

    statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))
    switch(method,
           "MC"={
               zzz <- .CI_MC_prim(x, statistic=statistic, B=10*B, slice.info=slice.info)
           },
           "SMC"={
               zzz <- .CI_MC_prim(x, statistic=statistic, B=B,    slice.info=slice.info)
               tot <- as.numeric(zzz[c("n.extreme","B")])
               if (slice.info){
                   slice <- zzz$slice
               }
               repeat{
                   if (tot[1]>L)
                       break
                   zzz<- .CI_MC_prim(x, statistic=statistic, B=B, slice.info=slice.info)
                   tot <- tot + as.numeric(zzz[c("n.extreme","B")])
                   if (slice.info){
                       slice[,"n.extreme"] <- slice[,"n.extreme"] + zzz$slice[,"n.extreme"]
                   }
               }
               
               zzz[c("p.value", "n.extreme", "B")] <- c(tot[1] / tot[2], tot)
               zzz["method"] <- "SMC"
               
               if (slice.info){
                   slice[,"p.value"] <- slice[,"n.extreme"] / tot[2]
                   zzz[["slice"]]    <- slice
               }
           })
    zzz
}


###
### CIP test; exact, based on monte carlo
###
.CI_MC_prim <- function(x, statistic="DEV", B=100, slice.info=FALSE){

    statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))
    
    .devFun2 <- function(obs, fit){ # Calculates deviance for independence model in r x c table
        ii  <- obs * fit > 0
        2 * sum(obs[ii] * log(obs[ii] / fit[ii]))
    }
    
    .X2Fun2 <- function(obs, fit){ # Calculates deviance for independence model in r x c table
        ii  <- obs * fit > 0
        a   <- (obs - fit)^2 / fit
        sum(a[ii])
    }

    .statFun <- if (statistic=="DEV") .devFun2 else .X2Fun2
    
    dn     <-  dim(x)
    v.idx  <-  seq_len(length(dn))
    v1R    <-  v.idx[-2]
    v2R    <-  v.idx[-1]
    dim12  <-  dim(x)[1:2]
    dim.R  <-  prod(dn[-(1:2)]) ## Careful when R is empty
    
    ## Marginal tables for (v1,R) and (v2,R) as matrices. Each row
    ## is a configuration of R
    t1R  <-  tableMargin(x, v1R)
    t2R  <-  tableMargin(x, v2R)
    t1R  <-  matrix(t1R, nrow=dim.R, byrow=TRUE)
    t2R  <-  matrix(t2R, nrow=dim.R, byrow=TRUE)
    xmat <-  matrix(x,   nrow=dim.R, byrow=TRUE)
    
    ## Find observed statistics for each slice
    tobs.slice <- vector("numeric", dim.R)
    for (ii in seq_len(dim.R)){
        r.sum    <- t1R[ii, ]
        c.sum    <- t2R[ii, ]
        expected <- outerPrim(r.sum, c.sum) / sum(r.sum)
        mm       <- xmat[ii, ]
        dim(mm)  <- dim12
        tobs.slice[ii] <- .statFun(mm, expected)
    }
    
    ## Find reference distribution for each slice
    tref.slice      <- matrix(NA, nrow=dim.R, ncol=B)
    n.extreme.slice <- vector("numeric", dim.R)
    for (ii in seq_len(nrow(t1R))){
        r.sum    <- t1R[ii,]
        c.sum    <- t2R[ii,]
        expected <- outerPrim(r.sum, c.sum) / sum(r.sum)
        zzz      <- r2dtable(B, r.sum, c.sum)
        for (kk in seq_len(B))
            tref.slice[ii, kk] <- .statFun(zzz[[kk]],expected)
        n.extreme.slice[ii] <- sum(tobs.slice[ii] < tref.slice[ii, ])
    }
    
    tref.total  <- colSums(tref.slice)
    tobs.total  <- sum(tobs.slice)
    n.extreme   <- sum(tobs.total < tref.total)
    p.value.slice <- n.extreme.slice / B
    p.value.total <- n.extreme / B
    
    if (slice.info){
        des   <- expand.grid(dimnames(x)[-(1:2)])
        slice <- cbind(data.frame(statistic=tobs.slice, n.extreme=n.extreme.slice,
                                  p.value=p.value.slice, df=NA), des)
    } else {
        slice=NULL
    }
    
    ans <- list(statistic=tobs.total, p.value=p.value.total, df=NA, statname=statistic,
                method="MC", varNames=names(dimnames(x)), n.extreme=n.extreme, B=B, slice=slice)
    class(ans) <- "citest"
    ans
}







