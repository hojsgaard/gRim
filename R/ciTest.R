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
#' @param set  A specification of the test to be made. The tests are of the form u
#'     and v are independent condionally on S where u and v are variables and S
#'     is a set of variables. See 'details' for details about specification of
#'     \code{set}.
#' @param \dots Additional arguments. 
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
#' @param set A specification of the test to be made. The tests are of the form
#'     u and v are independent condionally on S where u and v are variables and
#'     S is a set of variables. See 'details' for details about specification of
#'     \code{set}.
#' @param statistic The test statistic to be used, valid choices are
#'     \code{"DEV"} and \code{"F"}.
#' @param \dots Additional arguments
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


