###########################################################################
#'
#' @title Graphical Gaussian model
#' 
#' @description Specification of graphical Gaussian model. The 'c' in
#'     the name \code{cmod} refers to that it is a (graphical) model
#'     for 'c'ontinuous variables
#'
###########################################################################
#'
#' @details The independence model can be specified as \code{~.^1} and
#'     the saturated model as \code{~.^.}. The \code{marginal}
#'     argument can be used for specifying the independence or
#'     saturated models for only a subset of the variables.
#' @aliases extract_cmod_data
#' @param formula Model specification in one of the following forms:
#'     1) a right-hand sided formula, 2) as a list of generators.
#'     Notice that there are certain model specification shortcuts,
#'     see Section 'details' below.
#' 
#' @param data Data in one of the following forms: 1) A dataframe or
#'     2) a list with elements \code{cov} and \code{n.obs} (such as
#'     returned by the \code{cov.wt()} function.)
#'
#' @param marginal Should only a subset of the variables be used in
#'     connection with the model specification shortcuts.
#'
#' @param fit Should the model be fitted.
#'
#' @param details Control the amount of output; for debugging
#'     purposes.
#' @return An object of class \code{cModel} (a list)
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{dmod}}, \code{\link{mmod}},
#'     \code{\link{ggmfit}}
#' @keywords models
#' @examples
#' 
#' ## Graphical Gaussian model
#' data(carcass)
#' cm1 <- cmod(~ .^., data=carcass)
#' 
#' ## Stepwise selection based on BIC
#' cm2 <- backward(cm1, k=log(nrow(carcass)))
#' 
#' ## Stepwise selection with fixed edges
#' cm3 <- backward(cm1, k=log(nrow(carcass)),
#'  fixin=matrix(c("LeanMeat", "Meat11", "Meat12", "Meat13",
#'                 "LeanMeat", "Fat11", "Fat12", "Fat13"),
#'                  ncol=2))
#' 
#' @export cmod 
cmod <- function(formula, data, marginal=NULL, fit=TRUE, details=0) {

    cmod_data    <- extract_cmod_data(data)
    nms_in_data  <- colnames(cmod_data$S)

    model_spec   <- parse_gm_formula(formula, nms_in_data, marginal)
    nms_in_model <- model_spec$varNames
    glist        <- model_spec$glist

    ## Get varNames in the order matching to the data:
    nms_in_data <- nms_in_data[sort(match(nms_in_model, nms_in_data))]

    datainfo <- list(S     = cmod_data$S[nms_in_data, nms_in_data],
                     n.obs = cmod_data$n.obs,
                     data  = data)
    
    res <- list(modelinfo      = NULL,
                varNames       = nms_in_data,
                datainfo       = datainfo,
                fitinfo        = NULL,
                isFitted       = FALSE)
    
    upd   <- .cModel_finalize(glist, nms_in_data)  
    res$modelinfo <- upd  

    class(res) <- c("cModel", "iModel")
    if (fit)
        fit(res)
    else
        res
}

.cModel_finalize <- function(glist, varNames) {

    amat  <- ugList(glist, result="matrix")
    glist <- maxCliqueMAT(amat)[[1]]     ## FIXME: Rethink this
    isd   <- length(mcsMAT(amat)) > 0   
    
    list(glist       = glist,
         glistNUM    = .glistNUM(glist, varNames), 
         properties  = c(isg=TRUE, issd=isd))                
}


#' @export
fit.cModel <- function(object, engine=c("ncd", "covips", "conips"), start=NULL, ...) {
    dofit_cModel(object, engine=engine, start=start, ...)
}

## #' @export
dofit_cModel <- function(object, engine=c("ncd", "covips", "conips"), start=NULL, ...) {

    engine <- match.arg(engine)
    conips <- doBy::section_fun(wrapper_grips, list(method="conips"))
    covips <- doBy::section_fun(wrapper_grips, list(method="covips"))
    ncd    <- doBy::section_fun(wrapper_grips, list(method="ncd"))

    fitfun <- switch(engine,
                     ## "ggmfit"={ggmfit},
                     ## "ggmfitr"={ggmfitr},
                     "conips"={conips},
                     "covips"={covips},
                     "ncd"={ncd})
    
    ## Call C or R version of ips
    fit <- fitfun(object$datainfo$S,
                  n.obs = object$datainfo$n.obs,
                  glist = object$modelinfo$glist,
                  start = start, details=0,...)
    
    ## ideviance to independence model  
    idev  <-  fit$n.obs * (log(fit$detK) + sum(log(diag(fit$S))))  
    idim      <-  fit$nvar 
    sat.dim   <-  ((idim + 1) * idim) / 2
    dim.unadj <-  sat.dim - fit$df
    
    idf       <-  (dim.unadj - idim)
    logL.sat  <-  fit$logL + fit$dev / 2
    
    aic       <-  -2 * fit$logL + 2 * dim.unadj
    bic       <-  -2 * fit$logL + log(fit$n.obs) * dim.unadj
    
    dimension <- c(mod.dim=dim.unadj, sat.dim=sat.dim, i.dim=idim, df=fit$df, idf=idf)
    
    ans   <- list(dev=fit$dev, ideviance=idev, logL.sat=logL.sat,
                  aic=aic, bic=bic,
                  dimension=dimension)
    
    fit$S <- fit$n.obs <- fit$dev <- fit$df <- NULL
    ans  <- c(fit, ans)
    
    object$fitinfo  <- ans
    object$isFitted <- TRUE
    class(object)   <- c("cModel", "iModel")
    object
}

wrapper_grips <- function(S, n.obs, glist, start=NULL, 
                   eps=1e-12, iter=1000, details=0, method="ncd", ...){

    nms_in_data <- colnames(S)
    
    ## Numerical (indices) representation of glist
    glist.num <- lapply(glist, match, nms_in_data)
    
    fit   <- fit_ggm_grips(S, formula=glist.num, nobs=n.obs, method=method)
    detK  <- det(fit$K)
    dev   <- -n.obs * log(det(S %*% fit$K))            ## deviance to the saturated model  
    df    <-  sum(fit$K == 0) / 2
    
    out <- list(dev=dev, df=df, detK=detK, nvar=nrow(S), ## Pas på her
                S=S, n.obs=n.obs, logL=fit$details$logL,
                K=fit$K, iter=fit$details$iter)

    return(out)
}


#' @export
extract_cmod_data <- function(data.) {
    if (inherits(data., c("data.frame", "matrix"))) {
        data. <- cov.wt(data., method="ML")
    } else
        if (inherits(data., "list") &&
            identical(names(data.), c("cov", "center", "n.obs"))) {
            ## OK
        } else
            stop("Can not proceed...")
    
    names(data.)[1] <- "S"
    data.
}

