## Turn a formula into a list of generators (glist)
##
#' @title Parse graphical model formula
#'
#' @description Parse graphical model formula to internal representation
#'
#' @param formula A right hand sided formula or a list.
#' @param varnames Specification of the variables.
#' @param marginal Possible specification of marginal (a set of
#'     variables); useful in connection with model specification
#'     shortcuts.
#' @param interactions The maximum order of interactions allowed;
#'     useful in connection with model specification shortcuts.
parse_gm_formula <- function (formula, varnames, marginal=NULL, interactions=NULL)
{
    
    if (!inherits(formula, c("formula", "list"))) stop("Invalid formula specification")

    used.var <- if (length(marginal) > 0)  marginal  else  varnames 
    
    clformula <- class(formula) ##; cat("class(formula) :", clformula, "\n")
    
    switch(clformula,
           "formula"={
               pow <- .extract.power(formula)
               ##cat(sprintf("A formula is given; power=%d\n", pow))
               if (is.na(pow))
                   glist <- rhsFormula2list(formula) ##cat("A proper formula\n")
               else {
                   if (identical(pow, -1L)){
                       glist <- list(used.var)         ##cat("The saturated model\n")
                   } else {
                       if (identical(pow, 1L)){
                           glist <- as.list(used.var)  ##cat("The independence model\n")
                       } else {
                           pow   <- min(c(pow, length(used.var)))
                           glist <- gRbase::combnPrim(used.var, pow, simplify=FALSE)
                       }               
                   }
               }
           },
           "list"={
               glist <- formula
           })
    glist <- gRbase::removeRedundant(glist)    
    glist <- .check.glist(glist, used.var)  
    
    if (!is.null(interactions))
        glist <- .set.interactions(glist, interactions)
    
    value <- list(glist    = glist,
                  varNames = unique(unlist(glist)))
    value
}


.check.glist <- function(glist, used.var){


    if (any(is.na(pmatch(unlist(glist), used.var, 
                         duplicates.ok = TRUE)))) 
        stop("An invalid variable specification has been found\n")
    glist <- lapply(glist, function(x) {
        ii <- pmatch(x, used.var)
        used.var[ii]
    })
    
    modnames <- unique.default(unlist(glist))

    if (any(is.na(match(modnames, used.var))))
        stop("Variables in model not contained in the variable set. Perhaps a problem with 'marginal'?")
    
    glist
}


.set.interactions <- function(glist, interactions){
  zz <- lapply(glist, function(ss){
    if (length(ss) <= interactions){
        list(ss)
    } else {
        gRbase::combnPrim(ss, interactions, simplify=FALSE)
    }
  })
  
  gRbase::removeRedundant(unlist(zz, recursive=FALSE))
}



## .extract.power returns NA if the formula is "standard" or an
## integer if it has a "hat"
.extract.power <- function(form){

    if (length(form) == 3) stop("only rhs formula allowed")
    
    form.str <- deparse(form[[2]])

    has.hat <- length(grep("^\\.\\^", form.str)) > 0

    if (!has.hat) return (NA)

    rest <- gsub("\\.\\^", "", form.str)
    pp <- strsplit(rest, " ")[[1]][1]
    pow <- if (identical(pp, "."))
               -1L
           else
               as.integer(pp)
    
    if (identical(pow, 0L))
        pow <- 1L
    
    pow
}

