################################################################
##
## Test deletion of 'edge' from model 'object'
##
## If model is decomposable and edge is in one clique only, then
## degrees of freedom are adjusted for sparsity
##
################################################################

#' Test deletion of edge from an interaction model
#' 
#' Tests if an edge can be deleted from an interaction model.
#' 
#' If the model is decomposable and the edge is contained in one clique only
#' then the test is made in the marginal model given by that clique. In that
#' case, if the model is a log-linear model then degrees of freedom are
#' adjusted for sparsity
#' 
#' @aliases testdelete testdelete.iModel print.testdelete testdelete.mModel
#' @param object A model; an object of class \code{iModel}.
#' @param edge An edge in the model; either as a right-hand sided formula or as
#' a vector
#' @param k Penalty parameter used when calculating change in AIC
#' @param details The amount of details to be printed; 0 surpresses all
#' information
#' @param \dots Further arguments to be passed on to the underlying functions
#' for testing; that is to CItable and CImvn
#' @return A list.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{testadd}}
#' @keywords models htest
#' @examples
#' 
#' ## ## ## testdelete
#' ## ## ## 
#' 
#' ## ## Discrete model
#' ## ## 
#' data(reinis)
#' ## A decomposable model
#' ##
#' mf <- ~smoke:phys:mental + smoke:systol:mental
#' object <- dmod(mf, data=reinis)
#' 
#' testdelete(object, c("phys", "mental"))
#' testdelete(object, c("smoke", "mental"))
#' #testdelete(object, c("systol", "phys"))
#' 
#' ## A non-decomposable model
#' ##
#' mf <- ~smoke:phys + phys:mental + smoke:systol + systol:mental
#' object <- dmod(mf, data=reinis)
#' 
#' testdelete(object, c("phys", "mental"))
#' #testdelete(object, c("systol", "phys"))
#' #testdelete(object, c("smoke", "mental"))
#' 
#' ## ## Continuous model
#' ## ## 
#' data(math)
#' ## A decomposable model
#' ##
#' mf <- ~me:ve:al + me:al:an
#' object <- cmod(mf, data=math)
#' 
#' testdelete(object, c("ve", "al"))
#' testdelete(object, c("me", "al"))
#' 
#' @export testdelete
testdelete <- function(object, edge, k=2, details=1,...)
  UseMethod("testdelete")


testdelete.iModel <- function(object, edge, k=2, details=1, ...){

    model.type <- class(object)[1]
    ##cat(sprintf("testdelete.iModel model.type=%s\n", model.type))    

    edge <- rhsFormula2list(edge)[[1]]
    if (length(edge) !=2 )
        stop(paste("Not a valid edge: ", paste(edge, collapse=":"), " \n"))

    ## ----- START USING amat
    if (is.null((amat <- list(...)$amat)))
        amat <- .as_amat(getmi(object, "glist"))
        
    ## Is edge is in model? stop if not
    if (!subsetof(edge, colnames(amat)))
        stop(cat("variables:", edge, "not in model\n"))
    if (amat[edge[1], edge[2]] != 1)
        stop(cat("edge:", edge, "not in model\n"))
    
    ## Is model graphical?     ## FIXME: fails if model contains redundant elements..
    cliq    <- maxCliqueMAT(amat)$maxCliques
    isgraph <- length(cliq) == length(getmi(object, "glist"))

    ## Is model decomposable?
    isdecomp <- length(mcsMAT(amat)) > 0
    ## ----- STOP USING amat
    
    ## Is edge only in one clique in decomposable model?
    onlyinone <- FALSE
    if (isdecomp){
        idx   <- isin (cliq, edge, index=TRUE)
        onlyinone <- sum(idx) == 1
    }
    
    if (isdecomp && onlyinone && model.type %in% c("cModel", "dModel")){
        ## If edge is in one clique only, do test in marginal table
        hostcq <- cliq[idx == 1][[1]]
        set    <- c(edge, setdiff(hostcq, edge))
        ##cat(sprintf("CHK: edge: {%15s} hostcq: {%s}\n", toString(edge), toString(hostcq)))        
        ans <- switch(model.type,
                      "cModel"={ 
                          ciTest_mvn(list(cov=getmi(object, "S"),
                                          n.obs=getmi(object, "n")),
                                     set=set, ...)
                      },
                      "dModel"={
                          ciTest_table(getmi(object, "data"),
                                       set=set, ...)
                      })
        
        extra <- list(edge=edge, hostcq=hostcq, details=details, conmethod='data.based')
    } else {
        cat(sprintf("CHK: Make usual LR-test - edge = {%s}\n", toString(edge)))
        ob2   <- update(object, list(drop.edge=edge))
        ans   <- .comparemodels(object, ob2)
        extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
    }
    extra2 <- list(aic=ans$statistic - k * ans$df, k=k)
    ret    <- c(ans, extra, extra2)
    class(ret) <- "testdelete"
    ret
}

.comparemodels <- function(m1,m2) {
  devdiff <- 2 * (getmi(m2, "dev") - getmi(m1, "dev"))
  dfdiff  <- getmi(m2, "dimension")['df'] - getmi(m1, "dimension")['df']
  list('statistic'=devdiff, 'df'=dfdiff, 'p.value'=1 - pchisq(devdiff, dfdiff))
}

print.testdelete <- function(x,  ...){

  cat(sprintf("dev: %8.3f df:%3i p.value: %7.5f AIC(k=%3.1f): %6.1f edge: %s \n",
              x$statistic, x$df, x$p.value, x$k, x$aic, .toString(x$edge,':')))
  if (x$conmethod=="data.based"){
    if (x$details > 0){
      cat("host: ", x$hostcq, "\n")
      cat("Notice: Test performed in saturated marginal model\n")
    }
  } else {
    if (x$details > 0){
      cat("Notice: Test perfomed by comparing likelihood ratios\n")
    }
  }
  invisible(x)
}



