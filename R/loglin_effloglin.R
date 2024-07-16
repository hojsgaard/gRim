##########################################################################
#' @title Fitting Log-Linear Models by Message Passing
#' 
#' @description Fit log-linear models to multidimensional contingency
#'     tables by Iterative Proportional Fitting.
#'
#' @name loglin-effloglin
#' 
##########################################################################
#' 
#' @details The function differs from \code{loglin} in that 1) data
#'     can be given in the form of a list of sufficient marginals and
#'     2) the model is fitted only on the cliques of the triangulated
#'     interaction graph of the model. This means that the full table
#'     is not fitted, which means that \code{effloglin} is efficient
#'     (in terms of storage requirements). However \code{effloglin} is
#'     implemented entirely in R and is therefore slower than
#'     \code{loglin}. Argument names are chosen so as to match those
#'     of loglin()
#' @param table A contingency table
#' @param margin A generating class for a hierarchical log--linear model
#' @param fit If TRUE, the fitted values are returned.
#' @param eps Convergence limit; see 'details' below.
#' @param iter Maximum number of iterations allowed
#' @param print If TRUE, iteration details are printed.
#' @return A list.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{loglin}}
#'
#' @references Radim Jirousek and Stanislav Preucil (1995). On the
#'     effective implementation of the iterative proportional fitting
#'     procedure. Computational Statistics & Data Analysis Volume 19,
#'     Issue 2, February 1995, Pages 177-189
#' 
#' @keywords models
#' @examples
#' 
#' data(reinis)
#' glist <-list(c("smoke", "mental"), c("mental", "phys"),
#'              c("phys", "systol"), c("systol", "smoke"))
#' 
#' stab <- lapply(glist, function(gg) tabMarg(reinis, gg))
#' fv3 <- effloglin(stab, glist, print=FALSE)
#'
#' 
#' @export effloglin
effloglin <- function(table, margin, fit=FALSE, eps=0.01, iter=20, print=TRUE){

    ## margin is a list of generators
    xlogy <- function(x,y){
        i <- y > 0
        out <- rep(0, length(x))
        out[i] <- x[i] * log(y[i])
        out
    }
    
    amat <- ugList(margin, result="matrix")
    vn   <- colnames(amat)
    tri  <- triangulateMAT(amat)
    rip  <- ripMAT(tri)
    
    cliq     <- rip$cliques
    len.cliq <- length(cliq)
    
    ## get "host clique" for each generator
    ## FIXME: (effloglin) use general "get.host.clique" function.
    margin_host   <- rep(NA, length(margin))
    for (kk in 1:length(margin)){
        gg <- margin[[kk]]
        for (ii in seq_along(cliq)){
            zz <- match(gg, cliq[[ii]])
            if (!any(is.na(zz))){
                margin_host[kk] <- ii
                break
            }
        }
    }
    
    if (is.array(table)){
        Nobs   <- sum(table)
        suff_stat_list <-
            lapply(margin,
                   function(xx) {
                       tabMarg(table, xx)
                   })
    } else {
        Nobs   <- sum(table[[1]])
        suff_stat_list <- table
    }
    
    zzz       <- unlist(lapply(suff_stat_list, dimnames), recursive=FALSE)
    vl        <- zzz[unique.default(names(zzz))]

    cliq_pot  <- lapply(cliq, function(cq){
      gRain::cpt(cq, levels=vl[cq], values=1, normalize="all") 
    })
    
    ## Clique marginal probabilities
    cliq_prob  <- propagateLS(cliq_pot, rip, initialize=TRUE)    
        
    itcount  <- 1L
    logL     <- 0
    max.dif  <- vector("numeric", length(margin))
    repeat{
        for (ii in seq_along(margin)){        
            gg      <- margin[[ii]]
            cq      <- cliq[[margin_host[ii]]]
            cq.idx  <- margin_host[ii]
            cq_prob    <- cliq_prob[[cq.idx]]            
            st      <- suff_stat_list[[ii]]        ## Numerator in update
            tm      <- tabMarg(cq_prob, gg) * Nobs ## Denominator in update
            adjust2 <- tabDiv0(st, tm)             ## Adjustment
            max.dif[ii] <- max(abs(log(adjust2)))

            xly <- tabMult(st, log(adjust2))
            logL.adjust <- sum(xly, na.rm=TRUE)
            
            logL    <- logL + logL.adjust
            cat(glue("loop over margins. logL: {logL} \n\n"))
            
           cliq_pot[[cq.idx]] <- tabProd(cliq_pot[[cq.idx]], adjust2)
           cliq_prob          <- propagateLS(cliq_pot, rip, initialize=TRUE)
        }
      
        ## cat("iteration: ", itcount, "\n"); print(max.dif)
        if (print)
            cat("max deviation (obs-fitted):", max(max.dif), "\n")
        if ((max(max.dif) < eps) || (itcount >= iter))
            break()
        itcount <- itcount + 1L
    }
    
    vl    <- unlist(lapply(suff_stat_list, dimnames), recursive=FALSE)[vn]
    nlev  <- unlist(lapply(vl, length))
    gn    <- lapply(margin, match, vn)
    nparm <- .dim_loglin(gn, nlev)
    df    <- prod(nlev) - 1 - nparm
    
    ans <- list(potlist=cliq_pot, margin=margin, vn=vn, rip=rip, margin_host=margin_host,
                suff_stat_list=suff_stat_list, logL=logL, nparm=nparm, df=df)
    
### Create full joint:
    ##    print(cliq_prob)
    if (fit){
        pjoint <- cliq_prob[[1]]
        if (length(cliq_prob) > 1){
            for (ii in 2:length(cliq_prob)){
                sp <- rip$sep[[ii]]
                if (length(sp) > 0){              
                    mm <- tabMarg(cliq_prob[[ii]], sp)
                    cc <- tableOp(cliq_prob[[ii]], mm, "/")
                    pjoint <- tableOp(pjoint, cc, "*")
                } else {
                    pjoint <- tableOp(pjoint, cliq_prob[[ii]], "*")
                }
            }
        }
        pjoint <- tabPerm(pjoint, vn) * Nobs
        ans <- c(ans, list(fit=pjoint))
    }
    ## class(ans) <- "effloglin"
    return(ans)
}

