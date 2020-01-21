
#' @title General functions related to iModels
#'
#' @description General functions related to iModels
#'
#' @name imodel_general
#'
#' @param object,fit,x An \code{iModel} object.


#' @rdname imodel_general
logLik.iModel <- function(object, ...){
    val <- object$fitinfo$logL
    attr(val, "df") <- unname( object$fitinfo$dimension["mod.dim"] )
    attr(val, "nobs") <- sum(object$datainfo$data)
    class(val) <- "logLik"
    val
}


#' @rdname imodel_general
#' @param scale Unused (and irrelevant for these models)
#' @param k Weight of the degrees of freedom in the AIC formula
#' @param ... Currently unused.
#' 
extractAIC.iModel <- function(fit, scale, k = 2, ...){
    unname(c(fit$fitinfo$dimension["mod.dim"],
             -2*fit$fitinfo$logL + k*fit$fitinfo$dimension["mod.dim"]))
}

#' @rdname imodel_general
summary.iModel <- function(object, ...){
  glist <- getmi(object, "glist")
  isg   <- getmi(object, "isGraphical")
  isd   <- getmi(object, "isDecomposable")

  cq    <- getCliques(ugList(glist))# $maxCliques
  ans   <- structure(list(glist=glist, isGraphical=isg, isDecomposable=isd, cliques=cq),
                     class="iModelsummary")
  ans
}

#' @rdname imodel_general
print.iModelsummary <- function(x,...){
  cat(sprintf("is graphical=%s; is decomposable=%s\n", x$isGraphical, x$isDecomposable))
  cat("generators (glist):\n")
  str(x$glist, give.head=FALSE, comp.str=" ", no.list=TRUE)
  #cat("EXPERIMENTAL: components: ", names(x),"\n")
  invisible(x)
}

.extractFIT <- function(object,...){
  c(object[[1]], object$df)
}

.glist2formula <- function (f) {
  if (inherits(f, "formula"))
    return(f)
  ans <- try(as.formula(paste("~", paste(unlist(lapply(f, paste, collapse = "*")),
                                         collapse = "+")), .GlobalEnv),silent=TRUE)
  if (inherits(ans, "try-error"))
    stop("Unable to create formula from list. \nCould be due to white space, strange characters etc. in variable names\n")
  ans
}

#' importFrom stats formula terms

#' @export
#' @rdname imodel_general
formula.iModel <- function(x,...){
	#list2rhsFormula(x$glist)
  .glist2formula(x$glist)
}

#' @export
#' @rdname imodel_general
terms.iModel <- function(x, ...){
	x$glist
}


#' @rdname imodel_general
isGraphical.dModel <- function(x){
##    gRbase::isGraphical.default( terms(x) )
    isGraphical( terms(x) )
}

#' @rdname imodel_general
isDecomposable.dModel <- function(x){
    ##gRbase::isDecomposable.default( terms(x) )
    isDecomposable( terms(x) )
}

#' @rdname imodel_general
modelProperties <- function(object){
    UseMethod("modelProperties")
}

#' @rdname imodel_general                
modelProperties.dModel <- function(object){
    x <- terms( object )
    vn <- unique(unlist(x))
    amat <- glist2adjMAT(x, vn = vn)
    cliq <- maxCliqueMAT(amat)[[1]]
    isg <- all(unlist(lapply(cliq, function(cq) isin(x, cq))))
    isd <- if (isg) {
               length(mcsMAT(amat)) > 0
           }
           else FALSE
    
    c(isGraphical=isg, isDecomposable=isd)
}
