
#' logLik.iModel <- function(object,...)
#'   structure(object$fitinfo$logL, df=object$fitinfo$dimension["df"], class="logLik")

## Returns (df, AIC=-2logL + k df), so the objective is to mimimize this quantity
##

#' extractAIC.iModel <- function(fit, scale, k = 2, ...){
#'   unname(c(fit$fitinfo$dimension["df"],
#'            fit$fitinfo$dev - k*fit$fitinfo$dimension["df"]))
#' }


#' AIC.iModel <- function(object, ..., k=2){
#'     extractAIC(object)[2]
#' }


#' ## FIXME: need BIC for iModel
#' BIC.dModel <- function(object, ..., k=2){
#'     extractAIC(object, k=log(sum(object$datainfo$data)) )[2]
#' }


logLik.iModel <- function(object, ...){
    val <- object$fitinfo$logL
    attr(val, "df") <- unname( object$fitinfo$dimension["mod.dim"] )
    attr(val, "nobs") <- sum(object$datainfo$data)
    class(val) <- "logLik"
    val
}

extractAIC.iModel <- function(fit, scale, k = 2, ...){
    unname(c(fit$fitinfo$dimension["mod.dim"],
             -2*fit$fitinfo$logL + k*fit$fitinfo$dimension["mod.dim"]))
}

summary.iModel <- function(object, ...){
  glist <- object$glist

  isg   <- object$isGraphical
  isd   <- object$isDecomposable
  #cq    <- maxCliques(ugList(glist))$maxCliques
  cq    <- getCliques(ugList(glist))# $maxCliques
  ans   <- structure(list(glist=glist, isGraphical=isg, isDecomposable=isd, cliques=cq),
                     class="iModelsummary")
  ans
}

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

formula.iModel <- function(x,...){
	#list2rhsFormula(x$glist)
  .glist2formula(x$glist)
}

terms.iModel <- function(x,...){
	x$glist
}


isGraphical.dModel <- function(x){
    gRbase::isGraphical.default( terms(x) )
}

isDecomposable.dModel <- function(x){
    gRbase::isDecomposable.default( terms(x) )
}


modelProperties <- function(object){
    UseMethod("modelProperties")
}

modelProperties.dModel <- function(object){
    x <- terms( object )
    vn <- unique(unlist(x))
    amat <- glist2adjMAT(x, vn = vn)
    cliq <- maxCliqueMAT(amat)[[1]]
    isg <- all(unlist(lapply(cliq, function(sss) isin(x, sss))))
    isd <- if (isg) {
        length(mcsMAT(amat)) > 0
    }
    else FALSE
    c(isGraphical=isg, isDecomposable=isd)

}
