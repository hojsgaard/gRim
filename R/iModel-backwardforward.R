#########################################################################
##
## Backward and forward stepwise model selection for iModel's
##
#########################################################################

#' @title Stepwise model selection in (graphical) interaction models
#'
#' @description Stepwise model selection in (graphical) interaction models
#'
#' @name stepwise
#'
#' @param object An \code{iModel} model object
#' @param criterion Either \code{"aic"} or \code{"test"} (for significance test)
#' @param alpha Critical value for deeming an edge to be significant/
#'     insignificant. When \code{criterion="aic"}, \code{alpha} defaults to 0;
#'     when \code{criterion="test"}, \code{alpha} defaults to 0.05.
#' @param type Type of models to search. Either \code{"decomposable"} or
#'     \code{"unrestricted"}. If \code{type="decomposable"} and the initial
#'     model is decompsable, then the search is among decomposable models only.
#' @param search Either \code{'all'} (greedy) or \code{'headlong'} (search edges
#'     randomly; stop when an improvement has been found).
#' @param steps Maximum number of steps.
#' @param k Penalty term when \code{criterion="aic"}. Only k=2 gives genuine
#'     AIC.
#' @param fixin Matrix (p x 2) of edges. If those edges are in the model,
#'     they are not considered for removal.
#' @param fixout Matrix (p x 2) of edges. If those edges are not in the model,
#'     they are not considered for addition.
#' @param direction Direction for model search. Either \code{"backward"} or
#'     \code{"forward"}.
#' @param details Controls the level of printing on the screen.
#' @param trace For debugging only
#' @param ... Further arguments to be passed on to \code{testdelete} (for
#'     \code{testInEdges}) and \code{testadd} (for \code{testOutEdges}).
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{cmod}}, \code{\link{dmod}}, \code{\link{mmod}},
#'     \code{\link{testInEdges}}, \code{\link{testOutEdges}}
#' @keywords models
#' 
#' @examples
#'
#' data(reinis)
#' ## The saturated model
#' m1 <- dmod(~.^., data=reinis)
#' m2 <- stepwise(m1)
#' m2

stepwise.iModel <- function(object,
                          criterion="aic",
                          alpha=NULL,
                          type="decomposable",
                          search="all",
                          steps = 1000,
                          k = 2,
                          direction="backward",
                          fixin=NULL,
                          fixout=NULL,
                          details=0,
                          trace=2, ...)
{

    direction  <- match.arg(direction, c("backward", "forward", "both"))
    criterion  <- match.arg(criterion, c("aic",        "test"))
    search     <- match.arg(search,    c("headlong",   "all"))
    
    if (isGSD_glist(.glist(object))[2]){
        type <- match.arg(type,      c("decomposable", "unrestricted"))
    } else {
        type <- "unrestricted"
    }
    
    if(details>=1){
        cat("STEPWISE: ",
            "\n criterion:", criterion, if (criterion=="aic") paste("( k =", round(k,2),")"),
            "\n direction:", direction,
            "\n type     :", type,
            "\n search   :", search,
            "\n steps    :", steps, "\n")
    }
    
    if (direction == "backward"){
        ans <- backward(object, type=type, search=search,
                        criterion=criterion, alpha=alpha, steps=steps,
                        k=k, fixin=fixin, details=details, trace=trace)#$object
    } else {
        ans <- forward(object, type=type,  search=search,
                       criterion=criterion, alpha=alpha, steps=steps,
                       k=k, fixout=fixout, details=details, trace=trace)#$object
    }
    return(ans)
}

#' @rdname stepwise
backward <- function(object, criterion="aic", alpha=NULL, type="decomposable", search="all",
                     steps=1000,  k=2, fixin=NULL, details=1, trace=2, ...)
{
  type   <- match.arg(type,   c("decomposable", "unrestricted"))
  search <- match.arg(search, c("headlong",     "all"))
  ##details=10
  
  ## discrete is needed because checking for decomposability is special for mixed models.
  disc <- getmi(object, "disc.names") ## FIXME: Simpler than above
  vn   <- getmi(object, "varNames")

  fmat <- if (!is.null(fixin)) .as_amat(.do_handlefix(fixin), vn=vn)


  ## Here we use amat behind the scenes...
  isgsd  <- isGSD_glist(.glist(object), discrete=disc)
  if (!all(isgsd)) type <- "unrestricted"

  if (is.null(alpha))
    alpha <- if (criterion=="aic") 0 else 0.05
  
  switch(criterion,
       "aic" ={opt.op    <- which.min
               comp.op   <- `<`
               outstring <- "change.AIC"
               crit.str  <- "aic"},
       "test"={opt.op    <- which.max
               comp.op   <- `>`
               outstring <- "p.value"
               crit.str  <- "p.value"
               })

  t0 <- proc.time()

  testEdges <- switch(search,
                     "headlong" = {.testInEdges_headlong},
                     "all"      = {.testInEdges_all}) 
  itcount <- 1
  repeat{
      ##cat("Iteration", itcount, "\n")
      amat    <- .as_amat(.glist(object))
      edgeMAT <- getEdges(amat, type=type, ingraph=TRUE, discrete=disc) 

      if (!is.null(fmat))
          edgeMAT <- .update_edgeMAT(edgeMAT, fmat, vn)
      
      if (nrow(edgeMAT) == 0){
          if (details >= 1) cat(sprintf("No edges can be removed\n"))
          break
      } 
      
      testMAT <- testEdges(object, edgeMAT, comp.op=comp.op, crit.str=crit.str,
                           alpha=alpha, k=k, amat=amat, ...)
      statvec   <- testMAT[, crit.str]
      opt.idx   <- opt.op(statvec)

      if (details>=2) print(testMAT, row.names=FALSE, digits=4)
      ##str(list(statvec=statvec, opt.idx=opt.idx))
      if (comp.op( statvec[opt.idx], alpha)) {
          ## update object
          opt.edge <- as.character(testMAT[opt.idx, c("V1", "V2")])          
          object   <- update(object, list(drop.edge=opt.edge))
          if (details >= 1) cat(sprintf("  %s %9.4f Edge deleted: %s\n",
                                      outstring, statvec[opt.idx], .toString(opt.edge)))          
          if (itcount == steps) { break }
      } else break
      itcount <- itcount + 1
  }
  object
}


#' @rdname stepwise
forward <- function(object, criterion="aic", alpha=NULL, type="decomposable", search="all",
                     steps=1000,  k=2, fixout=NULL, details=1, trace=2, ...)
{
    type   <- match.arg(type,   c("decomposable", "unrestricted"))
    search <- match.arg(search, c("headlong",     "all"))
    
    ## discrete is needed because checking for decomposability is special for mixed models.
    disc <- getmi(object, "disc.names") ## FIXME: Simpler than above
    vn   <- getmi(object, "varNames")

    fmat <- if (!is.null(fixout)) .as_amat(.do_handlefix(fixout), vn=vn)
    
    isgsd  <- isGSD_glist(.glist(object), discrete=disc)    
    if (!all(isgsd))  type <- "unrestricted"
    
    if (is.null(alpha))
        alpha <- if (criterion=="aic") 0 else 0.05
    
    switch(criterion,
           "aic" ={opt.op    <- which.min
               comp.op   <- `<`
               outstring <- "change.AIC"
               crit.str  <- "aic"},
           "test"={opt.op    <- which.min
               comp.op   <- `<`
               outstring <- "p.value"
               crit.str  <- "p.value"
           })

    t0 <- proc.time()
    
    testEdges <- switch(search,
                      "headlong" ={.testOutEdges_headlong},
                      "all"      ={.testOutEdges_all})
    
    itcount <- 1
    repeat{
        ##cat("Iteration", itcount, "\n")
        amat    <- .as_amat(.glist(object))
        edgeMAT <- getEdges(amat, type=type, ingraph=FALSE, discrete=disc)  
        ##cat("edgeMAT:\n"); print(edgeMAT)
        
        if (!is.null(fmat))
            edgeMAT <- .update_edgeMAT(edgeMAT, fmat, vn)
        
        if (nrow(edgeMAT) == 0){
            if (details >= 1) cat(sprintf("No edges can be added\n"))
            break
        } 
        
        testMAT   <- testEdges(object, edgeMAT, comp.op=comp.op, crit.str=crit.str,
                               alpha=alpha, k=k, amat=amat, ...)
        
        statvec   <- testMAT[,crit.str]
        opt.idx   <- opt.op(statvec)
        
        if (details >= 2) print(testMAT, row.names=FALSE, digits=4)        
        if (comp.op( statvec[opt.idx], alpha)) {
            opt.edge    <- as.character(testMAT[opt.idx, c("V1","V2")])
            ## Update model 
            object  <- update(object, list(add.edge=opt.edge))            
            if (details >= 1) cat(sprintf("  %s %9.4f Edge added: %s\n",
                                        outstring, statvec[opt.idx], .toString(opt.edge)))
            if (itcount == steps) { break }
        } else break
        itcount <- itcount + 1
    }
    object
}




.amat_subtract <- function(amat1, amat2){
    if (is.null(amat2)) return(amat1)
    amat1 <- amat1 - amat2
    amat1[amat1 < 0L] <- 0L
    amat1
}

## Requires that first element (row) is smaller than second element (row)
.hash_edge <- function(em){
    z <- c(1L, 100000L)
    if (is.matrix(em) && nrow(em) == 2) colSums(em * z)
    else if (is.numeric(em) && length(em) == 2)sum(em * z)
}


.edge_matrix <- function(amat, duplicates=FALSE, names=FALSE, long=FALSE){
    if (!is.adjMAT(amat)) stop("amat not an adjacency matrix\n")
    m <- t.default(MAT2ftM_(amat)) ## Wide format
    if (!duplicates){
        i <- m[1, ] < m[2, ]
        m <- m[, i]
    }
    if (names) {
        d <- dim(m)
        m <- colnames(amat)[m]
        dim(m) <- d
    }
    if (long) t.default(m) else m
}

# obs: A named array
.iloglin <- function(obs, ...){
    vn <- names(dimnames(obs))
    fit <- lapply(vn, function(v) {ar_marg(obs, v)})
    n <- sum(obs)
    s <- lapply(fit, function(f) sum(f * log(f)))
    tt <- obs * log(obs)
    
    t0 <- sum(tt[!is.na(tt)])
    t1 <- sum(unlist(s))
    t2 <- length(fit) * n * log(n)
    dev <- 2 * (t0 - n * log(n) - t1 + t2)
    list(lrt=dev,
         df=prod(dim(obs)) - length(dim(obs)) - 1)
}


# x : generator liste
.is_gsd <- function(x, amat=NULL, discrete=NULL){
    if (is.null(amat)) amat <- .glist2amat(x)
    cliq <- maxCliqueMAT(amat)[[1]]
    ## Is graphical ?
    isg <- all(unlist(lapply(cliq, function(cq) isin(x, cq))))
    if (!isg) c(isg = FALSE, issd = FALSE)
    else
        c(isg = TRUE,
          issd = length(mcsmarkedMAT(amat, discrete = discrete)) > 0)
} 


.do_handlefix <- function(x){
    .handle_matrix <- function(x){
        if (ncol(x) == 2) rowmat2list(x)
        else if (nrow(x) == 2) colmat2list(x)
        else stop("Do not know what to do\n")        
    }
    
    .handle_list <- function(x){
        out <- lapply(x, function(l){
            if (is.vector(l)) list(l)
            else if (is.matrix(l)) .handle_matrix(l)
            else stop("Do not know what to do\n")        
        })
        unlist(out, recursive=FALSE)
    }
    
    if (is.null(x)) NULL
    else if (is.matrix(x)) .handle_matrix(x)
    else if (is.list(x)) .handle_list(x)
    else stop("Do not know what to do\n")
}

.update_edgeMAT <- function(edgeMAT, fmat, vn){
    am   <- .as_amat(matrix2list(edgeMAT), vn=vn)
    am   <- .amat_subtract(am, fmat)
    edgeMAT <- .edge_matrix(am, names=TRUE, long=TRUE)
    edgeMAT
}














##        





##setdiff_edge <- function(x, y){
##    ## x, y: p x 2 and q x 2 matrices - no checks
##    if (!(is.matrix(x) && ncol(x)==2)) stop("x is not p x 2 matrix\n")
##    if (!(is.matrix(y) && ncol(y)==2)) stop("y is not q x 2 matrix\n")
##    vn <- unique.default(c(x, y))
##    xn <- pairs2num(x, vn)
##    yn <- pairs2num(y, vn)    
##    d <-  setdiff(xn, yn)
##    i <-  match(d, xn)
##    x[i, ]
##}
##




##    .infoPrint(details, 2,
##               if (criterion=="aic"){
##                   cat(sprintf("FORWARD: type=%s search=%s, criterion=%s(%4.2f), alpha=%4.2f \n",
##                               type, search, criterion, k, alpha ))
##               } else {
##                   cat(sprintf("FORWARD: type=%s search=%s, criterion=%s, alpha=%4.2f \n",
##                               type, search, criterion, alpha ))
##               }
##               )
##    .infoPrint(details, 2,
##               cat(sprintf("Initial model: is graphical=%s is decomposable=%s\n",
##                           isgsd[1], isgsd[2])))
##




























