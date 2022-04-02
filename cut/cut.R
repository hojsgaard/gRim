



















































############################################################################







## #############################################################################


## testadd.iModel <- function(object, edge, k=2, details=1, ...){

##     model.type <- class(object)[1]
##     ##cat(sprintf("testadd.iModel model.type=%s\n", model.type))
    
##     edge <- rhsFormula2list(edge)[[1]]
##     if (length(edge)!=2)
##         stop(paste("Not a valid edge: ", paste(edge, collapse=":"), " \n"))
    
##     ## ----- START USING amat
##     if (is.null((amat <- list(...)$amat)))
##         amat <- .glist2adjMAT(object$glist)

##     ## Is edge is in model? stop if not
##     if (!subsetof(edge, colnames(amat)))
##         stop(cat("variables:", edge, "not in model\n"))

##     if (amat[edge[1],edge[2]]!=0)
##         stop(cat("edge:", edge, "already in model\n"))

##     ## Add edge to model FIXME: Fails if amat is sparse!
##     amat[edge[1], edge[2]] <- amat[edge[2], edge[1]] <- 1L

##     ## Is model graphical?
##     cliq <- maxCliqueMAT(amat)$maxCliques
##     isgraph <- length(cliq) == length(object$glist)
    
##     ## Is model decomposable?
##     isdecomp <- length(mcsMAT(amat)) > 0
##     ## ----- STOP USING amat
    
##     ## Is edge only in one clique in decomposable model?
##     onlyinone <- FALSE
##     if (isdecomp){
##         idx   <- isin (cliq, edge, index=TRUE)
##         onlyinone <- sum(idx) == 1
##     }

##   if (isdecomp && onlyinone && model.type %in% c("cModel","dModel")){
##       ## If edge is in one clique only, do test in marginal table
##       ##
##       hostcq <- cliq[idx==1][[1]]
      
##       set <- c(edge, setdiffPrim(hostcq, edge))
      
##       ans <- switch(model.type,
##                     "cModel"={ 
##                         ciTest_mvn(list(cov=getmi(object, "S"),
##                                         n.obs=getmi(object, "n")),
##                                    set=set, ...)
##                     },
##                     "dModel"={
##                         ciTest_table(getmi(object, "data"),
##                                      set=set, ...)
##                     })
            
##       extra <- list(edge=edge, hostcq=hostcq, details=details, conmethod='data.based')
##   } else {
##       ## Make usual LR-test
##       ##
##       ob2   <- update(object, list(add.edge=edge))
##       ans   <- .comparemodels(ob2,object)
##       extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
##   }
##     extra2 <- list(aic=-(ans$statistic - k * ans$df), k=k)
##     ret <- c(ans, extra, extra2)
##     class(ret) <- "testadd"
##     return(ret)
## }







##     ## cat("ans--------------\n"); print(ans)         
##     ## stop()
##     ## cat("ans2-------------\n");
##     ## ans2 <- .test_delete_edge(object, edge, details)
##     ## print(ans2)

##     ## ret <- .finalize_test(ans, k)
##     ## str(ret)
    
##     ## ans <- .test_delete_edge(object, edge, details)
##     ## oo <<- object
##     ## ee <<- edge
##     ## aa <<- ans



## .comparemodels <- function(m1, m2) {
##   devdiff <- 2 * (getmi(m2, "dev") - getmi(m1, "dev"))
##   dfdiff  <- getmi(m2, "dimension")['df'] - getmi(m1, "dimension")['df']
##   str(list(devdiff=devdiff, dfdiff=dfdiff))
##   list('statistic'=devdiff, 'df'=dfdiff, 'p.value'=1 - pchisq(devdiff, dfdiff))
## }




## testdelete.iModel <- function(object, edge, k=2, details=1, ...){

##     model.type <- class(object)[1]
##     ##cat(sprintf("testdelete.iModel model.type=%s\n", model.type))    

##     edge <- rhsFormula2list(edge)[[1]]
##     if (length(edge) !=2 )
##         stop(paste("Not a valid edge: ", paste(edge, collapse=":"), " \n"))

##     if (is.null((amat <- list(...)$amat)))
##         amat <- .as_amat(getmi(object, "glist"))
        
##     ## Is edge is in model? stop if not
##     if (!subsetof(edge, colnames(amat)))
##         stop(cat("variables:", edge, "not in model\n"))

##     if (amat[edge[1], edge[2]] != 1)
##         stop(cat("edge:", edge, "not in model\n"))
    
##     ## Is model graphical?     ## FIXME: fails if model contains redundant elements..
##     cliq    <- maxCliqueMAT(amat)$maxCliques
##     isgraph <- length(cliq) == length(getmi(object, "glist"))

##     ## Is model decomposable?
##     isdecomp <- length(mcsMAT(amat)) > 0
    
##     ## Is edge only in one clique in decomposable model?
##     onlyinone <- FALSE
##     if (isdecomp){
##         idx   <- isin (cliq, edge, index=TRUE)
##         onlyinone <- sum(idx) == 1
##     }
    
##     if (isdecomp && onlyinone && model.type %in% c("cModel", "dModel")){
##         ## If edge is in one clique only, do test in marginal table
##         hostcq <- cliq[idx == 1][[1]]
##         set    <- c(edge, setdiff(hostcq, edge))
##         ##cat(sprintf("CHK: edge: {%15s} hostcq: {%s}\n", toString(edge), toString(hostcq)))        
##         ans <- switch(model.type,
##                       "cModel"={ 
##                           ciTest_mvn(list(cov=getmi(object, "S"),
##                                           n.obs=getmi(object, "n")),
##                                      set=set, ...)
##                       },
##                       "dModel"={
##                           ciTest_table(getmi(object, "data"),
##                                        set=set, ...)
##                       })
        
##         extra <- list(edge=edge, hostcq=hostcq, details=details, conmethod='data.based')
##     } else {
##         cat(sprintf("CHK: Make usual LR-test - edge = {%s}\n", toString(edge)))
##         ob2   <- update(object, list(drop.edge=edge))
##         ans   <- .comparemodels(object, ob2)
##         extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
##     }
##     extra2 <- list(aic=ans$statistic - k * ans$df, k=k)
##     ret    <- c(ans, extra, extra2)
##     class(ret) <- "testdelete"
##     ret
## }




## ###### SEEMS UNUSED ######

## .preprocess <- function(e){
##     if( class(e)=="data.frame" ){
##         e <- as.matrix( e )
##     }
##     if (class(e)=="matrix"){
##         if (ncol(e) !=2 )
##             stop("matrix of edges must have two columns\n")
##         e <- rowmat2list( e )
##     }
##     if (!is.list(e))
##         e <- list(e)
##     e
## }

## addEdge_glist <- function(glist, e){
##     e <- .preprocess( e )
##     for (i in seq_along(e))
##         glist <- .add.edge_glist( glist, e[[i]] )
##     glist
## }

## dropEdge_glist <- function(glist, e){
##     e <- .preprocess( e )
##     for (i in seq_along(e))
##         glist <- .drop.edge_glist( glist, e[[i]] )
##     glist
## }

## addTerm_glist <- function(glist, e){
##     if (!is.list(e))
##         e <- list(e)
##     for (i in seq_along(e))
##         glist <- .add.term_glist( glist, e[[i]] )
##     glist
## }

## dropTerm_glist <- function(glist, e){
##     if (!is.list(e))
##         e <- list(e)

##     for (i in seq_along(e))
##         glist <- .drop.term_glist( glist, e[[i]] )
##     glist
## }

## ###### END ######














