################################################################################
## INITIAL VALUE
################################################################################


smallest_first_ordering <- function(adj_matrix) {
  # Calculate the degrees of each node (number of neighbors)
  degrees <- rowSums(adj_matrix)

  ## Initialize an empty vector to store the smallest-last ordering
  ordering <- numeric(length = ncol(adj_matrix))
  passive <- numeric(length = ncol(adj_matrix))

  for (i in 1:length(ordering)) {
    ## Find the node with the smallest degree
    degrees
    smallest_degree_node <- which(degrees == min(degrees))[1]
    smallest_degree_node

    ## Update the ordering and remove the node and its edges
    ordering[i] <- smallest_degree_node
    adj_matrix[smallest_degree_node,] <- 0
    adj_matrix[,smallest_degree_node] <- 0
    degrees <- rowSums(adj_matrix)
    degrees[smallest_degree_node] <- Inf  # Set the degree of the removed node to Inf
    passive[smallest_degree_node] <- Inf
    degrees <- degrees + passive
    ## str(list(passive=passive, degrees=degrees, ordering=ordering))

  }

  return(ordering)
}

small_first <- function(amat) {
  sfo <- smallest_first_ordering(amat)
  nms <- rownames(amat)
  sfo_inv <- match(1:nrow(amat), sfo)
  
  return(list(sfo=sfo, sfo_nms=nms[sfo],
              sfo_inv=sfo_inv,
              sfo_inv_nms=nms[sfo_inv]))
}

reorder <- function(S, amat) {
  tmp <- small_first(amat)
  S2 <- S[tmp$sfo, tmp$sfo]
  amat2 <- amat[tmp$sfo, tmp$sfo]
  sfo_inv <- match(1:nrow(amat), tmp$sfo)
  
  return(list(S2=S2,amat2=amat2, sfo=tmp$sfo, sfo_inv=sfo_inv))
}



################################################################################
################################################################################

r_ncd_ggm_ <- function(S, Elist, emat, nobs=NULL, K, iter, eps=1e-6, convcrit=1, print=FALSE, aux=NULL){
    t0 <- .get.time()

    ## cat(sprintf("iter: %d\n", iter))
    amat <- as_emat2amat(emat, nrow(S))
    parm <- new.env()

    parm$Sigma <- S
    parm$K <- NULL
    
    emat_c <- as_emat_complement(emat, nrow(S))    
    res1 <- outerloop1(parm, emat, emat_c, amat, eps, iter);
    ## cat("Sigma after outerloop1:\n"); print(parm$Sigma)
    
    res2 <- outerloop2(parm, emat, emat_c, amat, eps, iter)
    ## cat("Sigma and K after outerloop2 :\n"); print(parm$Sigma); print(parm$K)
    
    logL = ggm_logL_(parm$Sigma, K=parm$K, nobs=nobs)
    cert = certificate_ncd(parm$Sigma, parm$K, emat, emat_c)
    ## cat(sprintf("cert: %f\n", cert))

    out <- list(K=parm$K, Sigma=parm$Sigma,
                conv_check=unname(res2$conv_crit),
                logL = logL,
                cert = cert,
                iter = unname(res2$iter + res1$iter),
                time  = .get.diff.time(t0, units="millisecs"))
    return(out)
}


outerloop1 <- function(parm, emat, emat_c, amat, eps, maxit=1000){
    ## cat("outerloop1 - start\n")
    not_converged <- TRUE
    it1 <- 0

    if (!is.null(parm$K)){
        logLp <- ggm_logL_(parm$Sigma, parm$K, nobs=nobs)
    } else {
        logLp <- -99999
    }
    logL <- logLp

    Sigma_prev <- diag(diag(parm$Sigma))
    while (not_converged){
        innerloop1_update_Sigma2(parm, amat)
        it1  <- it1 + 1
        mad  <- mean_abs_diff_non_edge(parm$Sigma, Sigma_prev, emat_c)        
        Sigma_prev <- parm$Sigma
        conv_crit  <- mad

        ## cat(sprintf("it1: %d mad: %f\n", it1, mad))
        if ((it1 == maxit) || (conv_crit < eps)) { break() }
    }
    ## cat("outerloop1 done\n")
    list(iter=it1, mad=mad) ## FIXME mad should be conv_crit		
    ## parm (with updated Sigma) is the result 
}

outerloop2 <- function(parm, emat, emat_c, amat, eps, maxit=1000){

    smart_K <- TRUE
    
    is_invertible <- function(S){det(S) > 0}
    
    if (is_invertible(parm$Sigma)){
        parm$K <- solve_fun(parm$Sigma)        
    } else {
        stop("NCD algorithm failed")
    }

    ## cat("Sigma, K and KSigma (before updating)\n");
    ## print(parm$Sigma); print(parm$K); print(parm$K %*% parm$Sigma)
    
    it2       <- 0
    converged <- FALSE
    while (!converged){
        innerloop2_update_Sigma_K2(parm, amat, smart_K=smart_K)
        dif2 <- diff_fun(parm$Sigma, parm$K, emat_c)
        it2  <- it2 + 1
        conv_crit <- dif2
        ## cat(sprintf("it2: %d dif2: %f\n", it2, dif2))
        if ((it2 == maxit) || (conv_crit < eps)) { break() }
    }

    ## cat("outerloop2 done\n")
    out <- list(iter=it2, conv_crit=conv_crit)
    ## print(out)
    out
}





diff_fun <- function(S, K, emc){
    Kuv <- K[t(emc)]
    Suu <- S[cbind(emc[1,], emc[1,])]
    Svv <- S[cbind(emc[2,], emc[2,])]
    ## list(Kuv=Kuv, Suu=Suu, Svv=Svv, tr=sum(K*S), emc=emc, K=K) %>% print()
    d <- nrow(S) * (nrow(S) + 1) / 2
    sum( abs(Kuv) * sqrt(Suu * Svv) ) / d    
}

mean_abs_diff_non_edge <- function(Sigma1, Sigma2, emc){
    sum(abs(Sigma1 - Sigma2)) / (2 * ncol(emc))
}

certificate_ncd <- function(S, K, em, emc){
    Kuv <- K[t(emc)]

    K2 <- impose_zero(em, K)
    
    Suu <- S[cbind(emc[1,], emc[1,])]
    Svv <- S[cbind(emc[2,], emc[2,])]
    ## list(Kuv=Kuv, Suu=Suu, Svv=Svv, tr=sum(K*S), emc=emc, K=K) %>% print()
    sum(K2 * S) - nrow(S) + 2 * sum( abs(Kuv) * sqrt(Suu * Svv) ) 
    
}

## WITHOUT S
update_row_Sigma2 <- function(u, parm, amat, print=FALSE){

    ne_u2 <- amat[u, ] == 1 ## closure: u and its nbrs
    s_ru2 <- parm$Sigma[, u, drop = FALSE]    
    ## parm$Sigma ovenfor isf S
    
    if (all(!ne_u2)) {
        w <- rep(0, nrow(parm$Sigma))
    } else {
        beta <- rep(0, nrow(parm$Sigma)) 
        beta[ne_u2] <- solve_fun(parm$Sigma[ne_u2, ne_u2], s_ru2[ne_u2, ])
        w <- parm$Sigma %*% beta

        ## beta_star <- solve_fun(parm$Sigma[ne_u2, ne_u2], s_ru2[ne_u2, ])
        ## w_c <- parm$Sigma[bdu_c, bdu] %*% beta_star
    }

    tt <- parm$Sigma[u, u]
    parm$Sigma[, u] <- w ## s12 in Hastie notation
    parm$Sigma[u, ] <- w
    parm$Sigma[u, u] <- tt

    S2 <- parm$Sigma

    ## parm$Sigma[bdu_c, u] <- w_c ## s12 in Hastie notation
    ## parm$Sigma[u, bdu_c] <- w_c
}





## Without S
update_row_K2 <- function(u, parm, smart_K=TRUE, print=FALSE){

    if (smart_K){
        ## cat("smart update of K\n"); print(parm$K)
        BB <- parm$K[-u, u, drop=FALSE]
        ## cat("BB:\n"); print(t(BB))        
        CC_ <- parm$K[-u, -u]
        CC <- parm$K[-u, -u] - BB %*% t(BB) / as.numeric(parm$K[u, u])
        ## cat("CC_:\n"); print(CC_)
        ## cat("CC:\n"); print(CC)
        DD_ <- parm$Sigma[-u, u, drop=FALSE]
        DD <- CC %*% DD_ ##parm$Sigma[-u, u, drop=FALSE]
        ## cat("DD_:\n"); print(t(DD_))
        ## cat("DD:\n"); print(t(DD))
        k_uu_upd <- 1 / as.numeric(parm$Sigma[u, u] - parm$Sigma[u, -u, drop=FALSE] %*% DD) 

        ## NYT K
        parm$K[u, u]   <- k_uu_upd
        parm$K[-u, u]  <- -k_uu_upd * DD
        parm$K[u, -u]  <- t.default(-k_uu_upd * DD)        
        parm$K[-u, -u] <- CC + k_uu_upd * DD %*% t(DD)
        ## cat("K after smart update:\n"); print(parm$K)
    } else {
        ## cat("brute force update of K\n")
        parm$K <- solve_fun(parm$Sigma)        
    }    
}



## WITHOUT S
innerloop2_update_Sigma_K2 <- function(parm, amat, smart_K=TRUE, print=FALSE){
    ## cat("innerloop2_update_Sigma_K\n")
    for (u in 1:nrow(amat)){
        update_row_Sigma2(u, parm, amat, print=print)
        update_row_K2    (u, parm, smart_K=smart_K, print=print)
    }
}

## WITHOUT S
innerloop1_update_Sigma2 <- function(parm, amat){
    ## cat("innerloop1_update_Sigma\n")
    for (u in 1:nrow(amat)){
        update_row_Sigma2(u, parm, amat, print=FALSE)
    }
}



## outerloop2_s <- function(parm, emat, emat_c, amat, eps, maxit=1000){

##     is_invertible <- function(S){det(S) > 0}
    
##     if (is_invertible(parm$Sigma)){
##         parm$K <- solve_fun(parm$Sigma)        
##     } else {
##         stop("NCD algorithm failed")
##     }
    
##     it2       <- 0
##     dif2 <- diff_fun(parm$Sigma, parm$K, emat_c)
##     conv_crit <- dif2
    
##     ## cat("outerloop2 done\n")
##     out <- list(iter=it2, conv_crit=conv_crit)
##     ## print(out)
##     out
## }












