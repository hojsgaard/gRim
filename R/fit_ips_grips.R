
## --- ---------------------------------------------------
##
## ----- IMPLEMENTATION; ips + fips IPS for uggm
##
## --- ---------------------------------------------------

mean_abs_diff_on_emat <- function(Sigma, S, E){
    diff_on_emat <- function(Sigma, S, E){
        dif <- Sigma - S
        out <- c(diag(dif), dif[t(E)])
        out
    }

    mean(abs(diff_on_emat (Sigma, S, E)))
}


## #' @rdname fit-ggm
## #' @export
## edges can be either list or 2 x p matrix
r_covips_ggm_ <- function (S, Elist, emat, nobs, K, iter=1000L, eps=1e-6, convcrit=1, aux=list(), print=FALSE){

    t0 <- .get.time()
    inner <- get_inner_function(Elist)

    parm  <- get_init_parm(S, K)  
    nparm <- ncol(S) + ncol(emat)
    
    Scc_inv_list <- lapply(Elist, function(cc) {
        .sol(S[cc, cc, drop=FALSE])} )

    logLp <- ggm_logL_(S, K=parm$K, nobs=nobs)
    ## sprintf("logL (init): %f\n", logLp) %>% cat()
    
    it <- 0L
    repeat {
        it <- it + 1L        
        parm <- inner(Elist, S, parm, eps=eps, Scc_inv_list)

        switch(convcrit,
               "1" = {
                   conv_check = mean_abs_diff_on_emat(S, parm$Sigma, emat)
               },
               "2" = {
                   logL = ggm_logL_(S, K=parm$K, nobs=nobs)
                   conv_check <- abs((logL - logLp) / nparm)
                   logLp <- logL                   
               })        
        if (conv_check < eps || it == iter) break
    }

    if (convcrit==1)
        logL  = ggm_logL_(S, K=parm$K, nobs=nobs)

    parm <- c(parm,
              list(Sigma = solve_fun(parm$K),
                   iter  = it,
                   conv_check = conv_check,
                   time  = .get.diff.time(t0, units="millisecs"),
                   logL  = logL
                   ))        
    parm
}


get_inner_function <- function(edges){
    innerLoop_pair <- function(edges, S, parm, eps, Scc_inv_list){    
        for (j in 1:ncol(edges)){
            cc <- edges[, j]
            parm <- .covips_ggm_update_cc_parm(S, cc, parm, Scc_inv_list)
        }
        parm
    }
    
    innerLoop_gen <- function(edges, S, parm, eps, Scc_inv_list){    
        for (j in seq_along(edges)){
            cc <- edges[[j]]
            parm <- .covips_ggm_update_cc_parm(S, cc, parm, Scc_inv_list, j)
        }
        parm
    }

    if (is.matrix(edges))
        inner <- innerLoop_pair
    else if (is.list(edges))
        inner <- innerLoop_gen
    else stop("Can not find inner loop function")
    inner    
}


## edges is a list
r_conips_ggm_ <- function (S, Elist, emat, nobs, K, iter=1000, eps=1e-6, convcrit=1, aux=list(), print=FALSE) 
{
    my.complement <- function(cc, p){
        return(setdiff(1:p, cc))
    }

    nparm = ncol(S) + ncol(emat)
    
    ## t0 <- proc.time()
    t0 <- .get.time()
    if (!inherits(Elist, "list")) stop("'Elist' must be a list\n")

    if (length(Elist) == 1 && length(Elist[[1]]) == ncol(S)) ## Saturated model
    {
        out <- list(K    = solve_fun(S),
                    Sigma= S,
                    iter = 1,  ## Should be zero?
                    time = .get.diff.time(t0, units="millisecs"))
    }
    else if (length(Elist) == 0) ## Independence model
    {
        out <- list(K      = diag(1/diag(S)),
                    Sigma  = S,
                    iter   = 1, ## Should be zero?
                    time   =.get.diff.time(t0, units="millisecs"))
    }
    else
    {
        p  <- dim(S)[1]
        logLp = ggm_logL_(S, K=K, nobs=nobs)

        Scc_inv_list <-
            lapply(Elist, function(cc) {
                .sol(S[cc, cc, drop=FALSE])
            })
        
        Elist.complements <-
            lapply(Elist, my.complement, p)        

        it <- 0
        repeat {
            it <- it + 1
            for (j in 1:length(Elist)) {
                cc <- Elist[[j]]
                aa <- Elist.complements[[j]]
                D2 <- Scc_inv_list[[j]] +   ## solve(S[cc, cc, drop=FALSE]) +
                    K[cc, aa, drop=FALSE] %*%
                    solve_fun(K[aa, aa, drop=FALSE], K[aa, cc, drop=FALSE])                                
                K[cc, cc] <- D2
            }
            
            switch(convcrit,
                   "1" = {
                       conv_check = mean_abs_diff_on_emat(S, solve_fun(K), emat)
                   },
                   "2" = {
                       logL = ggm_logL_(S, K=K, nobs=nobs)
                       conv_check <- abs((logL - logLp) / nparm)
                       logLp <- logL                   
                   })
            
            if (conv_check < eps || it == iter) break
        }
        dimnames(K) <- dimnames(S)
        out <- list(K     = K,
                    Sigma = solve_fun(K),
                    conv_check = conv_check,
                    time  = .get.diff.time(t0, units="millisecs"),
                    iter  = it)
    }

    if (convcrit == 1)
        logL  = ggm_logL_(S, K=out$K, nobs=nobs)
   
    out$logL  = logL
    out
}

.covips_ggm_update_cc_parm <- function(S, cc, parm, Scc_inv_list, j)
{

    ## sprintf("parm-start:\n") %>% cat(); print(parm)
    
    Scc       <- S[cc, cc] ## FIXED
    Ktilde    <- Scc_inv_list[[j]]
    
    Kcc       <- parm$K[cc, cc]
    Sigmacc   <- parm$Sigma[cc, cc]
    
    Kstar  = .sol(Sigmacc);
    Laux   = Kcc - Kstar;
    
    Kupd   = Ktilde + Laux;
    Haux   = Kstar - Kstar %*% Scc %*% Kstar;
    
    parm$K[cc, cc] = Kupd;
    ## str(list(Kstar=Kstar, Haux=Haux, Scc=Scc))
    
    parm$Sigma     = parm$Sigma - parm$Sigma[, cc] %*% Haux %*% parm$Sigma[cc, ] 

    ## sprintf("parm-end:\n") %>% cat(); print(parm)
    parm
}







