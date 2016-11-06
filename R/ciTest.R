## ##################################################################
##
## Testing for condional independence CIP for various data types
## <x> : data
## <set>: NULL, a vector or a formula
##
## ##################################################################

#' @title Generic function for conditional independence test
#' 
#' @description Generic function for conditional independence test. Specializes
#'     to specific types of data.
#'
#' @name ciTest-general
#' 
#' @aliases ciTest ciTest.data.frame ciTest.table ciTest.list print.citest
#'     summary.citest
#' 
#' @param x An object for which a test for conditional independence is to be
#'     made. See 'details' for valid types of \code{x}.
#' @param set A specification of the test to be made. The tests are of the form
#'     u and v are independent condionally on S where u and v are variables and
#'     S is a set of variables. See 'details' for details about specification of
#'     \code{set}.
#' @param \dots Additional arguments to be passed on to other methods.
#' @return An object of class 'citest' (which is a list).
#'
#' @details \code{x} can be 1) a table, 2) a dataframe whose columns are
#'     numerics and factors or 3) a list with components \code{cov} and
#'     \code{n.obs}.
#'
#'     \code{set} can be 1) a vector or 2) a right-hand sided
#'     formula in which variables are separated by '+'. In either case, it is
#'     tested if the first two variables in the \code{set} are conditionally
#'     independent given the remaining variables in \code{set}.  (Notice an
#'     abuse of the '+' operator in the right-hand sided formula: The order of
#'     the variables does matter.)
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{ciTest.table}}, \code{\link{ciTest_table}},
#'     \code{\link{ciTest.data.frame}}, \code{\link{ciTest_df}},
#'     \code{\link{ciTest.list}}, \code{\link{ciTest_mvn}},
#'     \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#' 
#' ## contingency table:
#' data(reinis)
#' ## dataframe with only numeric variables:
#' data(carcass)
#' ## dataframe with numeric variables and factors:
#' data(milkcomp1)
#' 
#' ciTest(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
#' ciTest(reinis, set=~smo+phy+sys)
#' ciTest(milkcomp1, set=~tre+fat+pro)
#' 
#' 
#' @export ciTest
ciTest <- function(x, set=NULL, ...){
  UseMethod("ciTest")
}

ciTest.table <- function(x, set=NULL, ...){
  ciTest_table(x, set, ...)
}

ciTest.list <- function(x, set=NULL, ...){
  ciTest_mvn(x, set, ...)
}

ciTest.data.frame <- function(x, set=NULL, ...){
  ciTest_df(x, set, ...)
}

print.citest <- function(x, ...){
    if (length(x$varNames) > 2){
        cat("Testing", x$varNames[1], "_|_", x$varNames[2], "|",x$varNames[-(1:2)],"\n")
    } else {
        cat("Testing", x$varNames[1], "_|_", x$varNames[2], "\n")
    }
    cat(sprintf("Statistic (%s): %8.3f df: %s p-value: %6.4f method: %s\n",
                x$statname, x$statistic, x$df, x$p.value, x$method))

    if ( !is.null(x$slice) ){
        cat("Slice information:\n")
        print( x$slice, digits=4 )
    }

    invisible( x )
}

## FIXME: Do we need summary.citest??
summary.citest <- function(object,...){
    print( object )
    if ( !is.null(object$slice) ){
        cat("Slice information:\n")
        print( object$slice, digits=4 )
    }
    invisible( object )
}







## ########################################################
##
## CIP test in data.frame
## <x> : data.frame
##
## ########################################################

## ciTest_df <- function(x, set=NULL,...){
##   if (is.null(set)){
##     set <- names(x)
##   } else {
##     if (inherits(set,c("formula","character"))){
##       set <- unlist(rhsFormula2list(set))
##       set <- names(x)[pmatch(set, names(x))]
##     }
##   }

##   wdata       <- x[ ,set ]
##   varTypes    <- uniquePrim( unlist(lapply(wdata, class)) )

##   has.factor  <- "factor" %in% varTypes
##   has.numeric <- any(c("integer","numeric") %in% varTypes)
##   ##print(c(has.factor=has.factor, has.numeric=has.numeric))

##   if (has.factor & has.numeric){
##     .ciTest_df_internal(wdata, set,...)
##   } else {
##     if (has.factor){
##       ciTest_table(xtabs(~., data=wdata),set=set,...)
##     } else {
##       if (has.numeric){
##         ciTest_mvn(cov.wt(wdata,method="ML"), set=set,...)
##       } else {
##         stop("Strange error...\n")
##       }
##     }
##   }
## }



#' Test for conditional independence in a dataframe
#' 
#' Test for conditional independence in a dataframe.
#' 
#' \code{set} can be 1) a vector or 2) a right-hand sided formula in which
#' variables are separated by '+'. In either case, it is tested if the first
#' two variables in the \code{set} are conditionally independent given the
#' remaining variables in \code{set}.  (Notice an abuse of the '+' operator in
#' the right-hand sided formula: The order of the variables does matter.)
#' 
#' If \code{set} is \code{NULL} then it is tested whether the first two
#' variables are conditionally independent given the remaining variables.
#' 
#' If \code{set} consists only of factors then \code{x[,set]} is converted to a
#' contingency table and the test is made in this table using
#' \code{ciTest_table()}.
#' 
#' If \code{set} consists only of numeric values and integers then
#' \code{x[,set]} is converted to a list with components \code{cov} and
#' \code{n.obs} by calling \code{cov.wt(x[,set], method='ML')}. This list is
#' then passed on to \code{ciTest_mvn()} which makes the test.
#' 
#' @param x A dataframe.
#' @param set 
#' @param \dots 
#' @return An object of class 'citest' (which is a list).
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{ciTest}}
#' 
#' \code{\link{ciTest.table}} \code{\link{ciTest_table}}
#' 
#' % \code{\link{ciTest.data.frame}} % \code{\link{ciTest_df}}
#' 
#' \code{\link{ciTest.list}} \code{\link{ciTest_mvn}}
#' 
#' \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#' 
#' data(milkcomp1)
#' ciTest(milkcomp1, set=~tre+fat+pro)
#' ciTest_df(milkcomp1, set=~tre+fat+pro)
#' 
#' @export ciTest_df
ciTest_df <- function(x, set=NULL,...){
    if ( is.null( set ) ){
        set <- names( x )
    } else {
        if ( inherits( set, c("formula","character")) ){
            set <- unlist(rhsFormula2list(set), use.names = FALSE)
            set <- names( x )[ pmatch(set, names(x)) ]
        }
    }

    cat("set: "); print(set)

    wdata       <- x[ ,set ]
    varTypes    <- unique.default( unlist( lapply(wdata, class), use.names = FALSE) )

    has.factor  <- "factor" %in% varTypes
    has.numeric <- any(c("integer","numeric") %in% varTypes)
    switch.code <- as.character( 1*has.factor + 2*has.numeric )
    ##print(c(has.factor=has.factor, has.numeric=has.numeric))
    ##print(c(switch.code=switch.code))

    switch( switch.code,
           "0"={ ## F & F
               stop("Strange error...\n")
           },
           "1"={ ## T & F
               ciTest_table(xtabs(~., data=wdata),set=set, ...)
           },
           "2"={ ## F & T
               ciTest_mvn(cov.wt(wdata,method="ML"), set=set, ...)
           },
           "3"={ ## T & T
               .ciTest_df_internal(wdata, set, ...)
           }
           )
}




###
### If mixed data,  we test for deleting edge in model
###

### FIXME: set can be e.g. c(2,4,1) and this causes an error.
.ciTest_df_internal <- function(x, set=NULL,...){
    ##cat("CHK: ciTestmixed\n")
    if (is.numeric(set))
        set <- names(x)[ set ]

    obj <- mmod(list(set), data=x)

    ans <- testdelete(obj, set[1:2])
    ans <- ans[c(1,3,2)] ## FIXME: This is fragile
    ans$method   <- "CHISQ"
    ans$statname <- "DEV"
    ans$varNames <- set
    class(ans)   <- "citest"
    ans
}


## ########################################################
##
## CIP test in MVN-distribution
## <x>  : list(cov=, n.obs=)
##
## ########################################################



#' Test for conditional independence in the multivariate normal distribution
#' 
#' Test for conditional independence in the multivariate normal distribution.
#' 
#' \code{set} can be 1) a vector or 2) a right-hand sided formula in which
#' variables are separated by '+'. In either case, it is tested if the first
#' two variables in the \code{set} are conditionally independent given the
#' remaining variables in \code{set}.  (Notice an abuse of the '+' operator in
#' the right-hand sided formula: The order of the variables does matter.)
#' 
#' If \code{set} is \code{NULL} then it is tested whether the first two
#' variables are conditionally independent given the remaining variables.
#' 
#' \code{x} must be a list with components \code{cov} and \code{n.obs} such as
#' returned by calling \code{cov.wt( , method='ML')} on a dataframe.
#' 
#' @param x A list with elements \code{cov} and \code{n.obs} (such as returned
#'     from calling \code{cov.wt()} on a dataframe. See examples below.)
#' @param set
#' @param statistic
#' @param \dots
#' @return An object of class 'citest' (which is a list).
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{ciTest}} \code{\link{ciTest.table}},
#'     \code{\link{ciTest_table}}, \code{\link{ciTest.data.frame}},
#'     \code{\link{ciTest_df}}, \code{\link{ciTest.list}},
#'     \code{\link{ciTest_mvn}}, \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#' 
#' data(carcass)
#' ciTest(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
#' ciTest_mvn(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
#' 
#' @export ciTest_mvn
ciTest_mvn <- function(x,set=NULL, statistic="DEV",...){
  if(any(is.na(match(c("cov","n.obs"), names(x))))){
    stop("Expecting a list with components 'cov' and 'n.obs'\n")
  }

  if (is.null(set)){
    set   <- colnames( x$cov )
    x$cov <- x$cov[ set, set ]
  } else {
    if (inherits(set,c("formula","character"))){
      set <- unlist(rhsFormula2list(set))
      set <- colnames(x$cov)[pmatch(set, colnames(x$cov))]
      x$cov <- x$cov[set,set]
    }
  }
  .ciTest_mvn_internal(x, statistic=statistic, ...) ## Should <set> go here...
}

### This is the workhorse.
###
.ciTest_mvn_internal <- function(x, statistic="DEV", ...){

  statistic <- match.arg(toupper(statistic), c("DEV","F"))
  if (statistic=="DEV")
    method <- "CHISQ"
  else
    method <- "F"

  S     <- x$cov
  n.obs <- x$n.obs

  vn <- colnames(S)
  K  <- length(vn)

  R   <- vn[-(1:2)]
  v1R <- c(vn[1],R)
  v2R <- c(vn[2],R)

  v1R.idx <- match(v1R, vn)
  v2R.idx <- match(v2R, vn)
  R.idx   <- match(R,   vn)

  d <- n.obs * (log(det(S[v1R.idx, v1R.idx, drop=FALSE])) +
                log(det(S[v2R.idx, v2R.idx, drop=FALSE])) -
                log(det(S[R.idx, R.idx, drop=FALSE])) - log(det(S))
                )

  num.df <- 1
  switch(statistic,
         "DEV"={
           tobs     <- d
           denom.df <- NULL
           p        <- 1-pchisq(tobs, df=num.df)
         },
         "F"={
           tobs     <- (exp(d/n.obs)-1)*(n.obs-K)
           denom.df <- n.obs-K
           p        <- 1-pf(tobs, num.df, denom.df)
         })

  ans <- list(statistic=tobs, p.value=p, df=num.df, denom.df=denom.df,
              statname=statistic, method=method, varNames=vn)
  class(ans) <- "citest"
  ans
}

## ########################################################
##
## CIP test in a table
## <x>  : table
##
## ########################################################



#' Test for conditional independence in a contingency table
#' 
#' Test for conditional independence in a contingency table
#' 
#' \code{set} can be 1) a vector or 2) a right-hand sided formula in which
#' variables are separated by '+'. In either case, it is tested if the first
#' two variables in the \code{set} are conditionally independent given the
#' remaining variables in \code{set}.  (Notice an abuse of the '+' operator in
#' the right-hand sided formula: The order of the variables does matter.)
#' 
#' If \code{set} is \code{NULL} then it is tested whether the first two
#' variables are conditionally independent given the remaining variables.
#' 
#' @param x A contingency table.
#' @param set
#' @param statistic
#' @param method
#' @param adjust.df
#' @param slice.info
#' @param L
#' @param B
#' @param ...  Additional arguments.
#' @return An object of class 'citest' (which is a list).
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{ciTest}}, \code{\link{ciTest.data.frame}},
#'     \code{\link{ciTest_df}}, \code{\link{ciTest.list}},
#'     \code{\link{ciTest_mvn}}, \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#' 
#' data(reinis)
#' ciTest(reinis, set=~smo+phy+sys)
#' ciTest_table(reinis, set=~smo+phy+sys)
#' 
#' @export ciTest_table
ciTest_table <- function(x, set=NULL, statistic="dev", method="chisq", adjust.df=TRUE, slice.info=TRUE, L=20, B=200, ...){

  statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))
  method    <- match.arg(toupper(method),    c("CHISQ", "MC", "SMC"))

  if (is.null(set)){
    set <- names(dimnames(x))
  } else {
    if (inherits(set,"numeric")){
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

  t.v1R <- tableMargin(x, c(v1,R))

  t.v2R <- tableMargin(x, c(v2,R))

  ## Fit table
  if (length(R)){

    t.R   <- tableMargin(x, R)
    fit.table <- tablePerm(tableOp(tableOp(t.v1R, t.v2R), t.R, "/"), vn)

  } else {
    fit.table <- tablePerm(tableOp(t.v1R, t.v2R), vn)/sum(x)
  }

  ## Evaluate test statistic
  if (statistic=="DEV"){             ## Deviance
    tobs  <- 2* x * log(x/fit.table)
  } else {                           ## Pearson X2
    tobs <- (x-fit.table)^2/fit.table
  }
  tobs[!is.finite(tobs)] <- 0
  tobsGlobal <- sum(tobs)



  ## Calculate df with or without adjustment for sparsity
  if (adjust.df){
    t.v1Rmat    <- matrix(t.v1R, nrow=dim.R,byrow=TRUE)
    t.v2Rmat    <- matrix(t.v2R, nrow=dim.R,byrow=TRUE)
    zzz <- (t.v1Rmat>0)*1
    if (!is.null(dim(zzz)))
      dim.v1.adj  <- rowSumsPrim(zzz)
    else
      dim.v1.adj  <- sum(zzz)

    zzz <- (t.v2Rmat>0)*1
    if (!is.null(dim(zzz)))
      dim.v2.adj  <- rowSumsPrim(zzz)
    else
      dim.v2.adj  <- sum(zzz)

    d1        <- dim.v1.adj-1
    d1[d1<0]  <- 0
    d2        <- dim.v2.adj-1
    d2[d2<0]  <- 0
    dofSlice  <- d1 * d2
  } else {
    dofSlice <- rep.int((dim.v1-1)*(dim.v2-1), dim.R)
  }
  dofGlobal <- sum(dofSlice)
  pGlobal   <- 1-pchisq(tobsGlobal, dofGlobal)

  tobsSlice <- rowSumsPrim(matrix(tobs, nrow=dim.R, byrow=TRUE))
  pSlice    <- 1-pchisq(tobsSlice, df=dofSlice)

  if (length(R) && slice.info){
    sliceInfo <- list(statistic=tobsSlice, p.value=pSlice, df=dofSlice)
    des   <- expand.grid(dimnames(x)[-(1:2)])
    slice <- cbind(as.data.frame(sliceInfo[1:3]), des)
  } else {
    slice <- NULL
  }
  ans <- list(statistic=tobsGlobal, p.value=pGlobal, df=dofGlobal, statname=statistic, method="CHISQ",
              adjust.df=adjust.df, varNames=names(dimnames(x)),slice=slice)
  class(ans) <- "citest"
  ans
}


###
### CIP test; exact, based on sequential monte carlo
###

.CI_SMC_prim <- function(x, statistic="DEV", method="SMC", L=50, B=200, slice.info=FALSE){

  statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))

  switch(method,
         "MC"={
           zzz <- .CI_MC_prim(x, statistic=statistic,B=10*B, slice.info=slice.info)
         },
         "SMC"={
           zzz <- .CI_MC_prim(x, statistic=statistic, B=B, slice.info=slice.info)
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

           zzz[c("p.value","n.extreme","B")] <- c(tot[1]/tot[2], tot)
           zzz["method"] <- "SMC"

           if (slice.info){
             slice[,"p.value"] <- slice[,"n.extreme"]/tot[2]
             zzz[["slice"]]    <- slice
           }
         })
  return(zzz)
}


###
### CIP test; exact, based on monte carlo
###
.CI_MC_prim <- function(x, statistic="DEV", B=100, slice.info=FALSE){

  statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))

##   .r2dtable <- function(n,r,c){
##     #.Call("R_r2dtable", as.integer(n), as.integer(r), as.integer(c), PACKAGE = "base")
##     r2dtable(n,r,c)
##   }

  .devFun2 <- function(mm,expected){ # Calculates deviance for independence model in r x c table
    ii  <- mm * expected > 0
    2*sum(mm[ii]*log(mm[ii]/expected[ii]))
  }

  .X2Fun2 <- function(mm,expected){ # Calculates deviance for independence model in r x c table
    ii  <- mm * expected > 0
    a   <- (mm-expected)^2/expected
    sum(a[ii])
  }

  if (statistic=="DEV")
    .statFun <- .devFun2
  else
    .statFun <- .X2Fun2

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
     r.sum    <- t1R[ii,]
     c.sum    <- t2R[ii,]
     expected <- outerPrim(r.sum,c.sum)/sum(r.sum)
     mm       <- xmat[ii,]
     dim(mm)  <- dim12
     tobs.slice[ii] <- .statFun(mm, expected)
  }

  ## Find reference distribution for each slice
  tref.slice      <- matrix(NA, nrow=dim.R, ncol=B)
  n.extreme.slice <- vector("numeric", dim.R)
  for (ii in seq_len(nrow(t1R))){
    r.sum    <- t1R[ii,]
    c.sum    <- t2R[ii,]
    expected <- outerPrim(r.sum,c.sum)/sum(r.sum)
    zzz      <- r2dtable(B, r.sum, c.sum)
    for (kk in seq_len(B))
      tref.slice[ii,kk] <- .statFun(zzz[[kk]],expected)
    n.extreme.slice[ii] <- sum(tobs.slice[ii] < tref.slice[ii,])
  }

  tref.total  <- colSumsPrim(tref.slice)
  tobs.total  <- sum(tobs.slice)
  n.extreme   <- sum(tobs.total<tref.total)
  p.value.slice <- n.extreme.slice/B
  p.value.total <- n.extreme/B

  if (slice.info){
    des   <- expand.grid(dimnames(x)[-(1:2)])
    slice <- cbind(data.frame(statistic=tobs.slice, n.extreme=n.extreme.slice, p.value=p.value.slice, df=NA),des)
  } else {
    slice=NULL
  }

  ans <- list(statistic=tobs.total, p.value=p.value.total, df=NA, statname=statistic,
                       method="MC", varNames=names(dimnames(x)), n.extreme=n.extreme,B=B,slice=slice)
  class(ans) <- "citest"

  return(ans)
}







 
