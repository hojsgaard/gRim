#' @title Edge matrix operations
#' @description Edge matrix operations needed for ips algorithms
#' @name emat_operations
#' @note These functions may well be removed from the package in future relases
#' @details An emat with p edges is represented by a 2 x p matrix.
#'
#' @param emat,emat1,emat2 Edge matrix (a 2 x p matrix)
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @examples
#' emat1 <- model_saturated(3:4, type="emat")
#' emat2 <- model_saturated(1:4, type="emat")
#' emat_complement(emat1, emat2)
#' emat3 <- model_saturated(2:4, type="emat")
#' emat_compare(emat1, emat3)

#' @rdname emat_operations
#' @export
emat_compare <- function(emat1, emat2){
    ## emat1, emat2 : 2 x p matrices
    ## common: Edges in emat1 AND emat2
    ## emat1only: Edges in emat1 that are not in emat2
    ## emat2only: Edges in emat2 that are not in emat1    
    emat1 <- colmat2list(emat1)
    emat2 <- colmat2list(emat2)

    common=emat1[(sapply(emat1, is_inset, emat2))]
    emat1only=emat1[!sapply(emat1, is_inset, emat2)]
    emat2only=emat2[!sapply(emat2, is_inset, emat1)]

    common=do.call(cbind, common)
    emat1only=do.call(cbind, emat1only)
    emat2only=do.call(cbind, emat2only)
    out  <- list(common=common, emat1only=emat1only, emat2only=emat2only)
    out
}

#' @rdname emat_operations
#' @export
emat_complement <- function(emat1, emat2){
    ## emat1, emat2 : 2 x p matrices
    ## Find complement of emat1 in emat2
    gl1 <- colmat2list(emat1)
    gl2 <- colmat2list(emat2)
    emat2[, !sapply(gl2, is_inset, gl1)]
}

#' @rdname emat_operations
#' @export
emat_sort <- function(emat1){
    id <- emat1[1,] > emat1[2,]
    emat1[1:2, id] <- emat1[2:1, id]
    o <- order(emat1[1,])
    emat1[, o]
}

#' @rdname emat_operations
#' @export
order_rows <- function(emat){
    b <- emat[1,] > emat[2,]
    emat[1:2, b] <- emat[2:1, b]
    emat
}

## #################################################################

#' @title Generate various grapical models
#' @description Models are represented in various forms
#' @name generate_models
#' @param index A vector of integers
#' @param dim A vector with dimensions
#' @param prob Probability of any edge being present.
#' @param nms Names of variables.
#' @param type Output type.
#' 



emat_saturated_model <- function(index){
    t(names2pairs(index, result="matrix"))
}

#' @rdname generate_models
#' @export
model_saturated <- function(index, type="emat", nms=NULL){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_saturated_model(index)
    out <- switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )
    if (!is.null(nms))
        add_nms(out, nms, type)
    else
        out

}




## #################################################################

emat_random_tree <- function(index, prob=0){

    emat_random_tree0 <- function(index){
        ## Based on https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence
        N <- max(index)
        emat <- matrix(0, nrow=2, ncol=N-1)
        degree <- rep(1, N)
        S <- sample(N, N-2, replace=TRUE)
        for (i in S){
            degree[i] <- degree[i] + 1
        }
        
        k <- 1
        for (i in S){
            for (j in seq_len(N)){
                ## cat(sprintf("i %d j %d\n", i, j))
                if (degree[j] == 1){
                    edge <- c(j, i)
                    ## cat("edge: \n"); print(edge)
                    emat[, k] <- edge
                    k <- k + 1
                    degree[i] <- degree[i] - 1
                    degree[j] <- degree[j] - 1
                    break
                }
            }    
        }
        emat[, k] <- which(degree == 1)    
        emat
    }

    if ((prob < 0) || (prob > 1)) stop("prob must be in [0, 1]")
    if (prob < 1e-6)
        return(emat_random_tree0(index))

    emat1 <- emat_random_tree(index)
    emat2 <- emat_random_model(index, prob)
    ee   <- cbind(order_rows(emat1), order_rows(emat2))
    emat <- t.default(unique(t.default(ee)))
    emat
}

#' @rdname generate_models
#' @export
model_random_tree <- function(index, prob=0, type="emat", nms=NULL){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_random_tree(index, prob)
    out <- switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )

    if (!is.null(nms))
        add_nms(out, nms, type)
    else
        out
}


## #################################################################

emat_rectangular_grid <- function(dim){

    if (length(dim)==1)
        dim <- rep(dim, 2)

    nr.nc <- dim
    nr <- nr.nc[1]
    nc <- nr.nc[2]
    
    ## first row of edges
    rr <- rbind(
        1:(nc-1),
        2:nc
    )
    
    ## first column of edges
    bb <- 0:(nr-1)*nc + 1
    cc <- rbind(
        bb[1:(nr-1)],
        bb[2:nr]
    )
        
    ## all rows of edges
    roe <- lapply((0:(nr-1))*nc, function(r) rr + r)  |>  do.call(cbind, args=_)
    ## all columns of edges
    coe <- lapply((0:(nc-1)), function(c) cc + c)  |>  do.call(cbind, args=_)
    
    out <- cbind(roe, coe)
    out
}

#' @rdname generate_models
#' @export
model_rectangular_grid <- function(dim, type="emat", nms=NULL){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_rectangular_grid(dim)
    out <- switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, prod(dim))},
           "amat" ={as_emat2amat(em, prod(dim))}
           )

    if (!is.null(nms))
        add_nms(out, nms, type)
    else
        out

}

## #################################################################

emat_line_model <- function(index){
  unname(rbind(index[-length(index)], index[-1]))
}

#' @rdname generate_models
#' @export
model_line <- function(index, type="emat", nms=NULL){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_line_model(index)
    out <- switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )

    if (!is.null(nms))
        add_nms(out, nms, type)
    else
        out
}


## #################################################################

emat_star_model <- function(index){
    unname(rbind(index[1], index[-1]))
}

#' @rdname generate_models
#' @export
model_star <- function(index, type="emat", nms=NULL){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_star_model(index)
    out <- switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )

    if (!is.null(nms))
        add_nms(out, nms, type)
    else
        out
}


## #################################################################

emat_loop_model <- function(index, prob=0){

    emat_loop_model0 <- function(index){
        if (length(index) == 2) matrix(index)
        else unname(rbind(index, c(index[-1], index[1])))
    }
    
    if ((prob < 0) || (prob > 1)) stop("prob must be in [0, 1]")
    if (prob < 1e-6)
        return(emat_loop_model0(index))
    
    emat1 <- emat_loop_model0(index)
    emat2 <- emat_random_model(index, prob)
    ee   <- cbind(order_rows(emat1), order_rows(emat2))
    emat <- t.default(unique(t.default(ee)))
    emat
}

#' @rdname generate_models
#' @export
model_loop <- function(index, prob=0, type="emat",nms=NULL){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_loop_model(index, prob)
    out <- switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )

    if (!is.null(nms))
        add_nms(out, nms, type)
    else
        out

}


## #################################################################

emat_random_model <- function(index, prob=0.1){
    N <- max(index)
    g     <- igraph::sample_gnp(N, prob)
    emat  <- t(igraph::get.edgelist(g))
    emat    
}

#' @rdname generate_models 
#' @export
model_random <- function(index, prob=0.1, type="emat",nms=NULL){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_random_model(index, prob)
    out <- switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )

    if (!is.null(nms))
        add_nms(out, nms, type)
    else
        out
}




add_nms <- function(x, nms, type){    
    if (length(x) == 0) return(x)
    switch(type,
           "emat" ={if (ncol(x)==0){
                        x                                                
                    } else {
                        x2 <- nms[x]
                        dim(x2) <- dim(x)
                        x2
                    }},
           
                        
                        
                        
           "amat"={
               rownames(x)<-colnames(x)<-nms[1:nrow(x)]
               x
           },
           "glist"=,"cliq"={
               w <- lapply(x, function(z) {
                   nms[z]
               })
               w
           })
}    




#' @title Coerce models to different representations
#' @description Coerce models to different representations
#' @name coerce_models
#'
#' @param emat Edge matrix (2 x p)
#' @param amat Adjacency matrix
#' @param glist Generator list
#' @param elist Edge list (list of pairs)
#' @param nvar Number of variables
#' @param K Concentration matrix
#' @param d Number of columns in output.
#' @param eps Small number
#'
#' glist <- list(c(1,2,3),c(2,3,4),c(5,6))
#' em <- as_glist2emat(glist)
#' am <- as_emat2amat(em, d=6)
#' ig <- as_emat2igraph(em)
#' el <- as_emat2elist(em)
#' igraph::max_cliques(ig)
#' as_emat2cq(em, 6)
#' as_emat_complement(em, 6)


# Convert edges to cliques
#' @export
#' @rdname coerce_models
as_emat2cq <- function(emat, nvar=NULL){

    gg <- graph_from_edgelist(t(emat), directed=FALSE)
    gRbase::maxClique(gg)$maxCliques

    am <- as_adjacency_matrix(gg, sparse = FALSE)

    
    if (ncol(emat) == 0) return(list())
    if (is.null(nvar)){
        nvar <- max(emat)
    }
    am <- matrix(0, nrow=nvar, ncol=nvar)    
    am[t(emat)] <- 1  
    am <- am + t.default(am)

    cq <- as(am, "igraph") |> as("dgCMatrix") |> gRbase::maxCliqueMAT() 
    cq <- lapply(cq$maxCliques, "as.numeric")
    cq
}


#' @export
#' @rdname coerce_models
as_emat_complement <- function(emat, nvar){ # Those edges NOT in emat
  am <- matrix(0, nrow=nvar, ncol=nvar)
  am[t(emat)] <- 1
  am[t(emat[2:1, ,drop=FALSE])] <- 1
  am <- (am == 0) * lower.tri(am)
  t.default(which(am==1, arr.ind=TRUE))
}


#' @export
#' @rdname coerce_models
as_emat2amat <-function(emat, d){
    vn <- sort(unique(c(emat)))
    M <- matrix(0, nrow=d, ncol=d)
    if (ncol(emat)>0){
        for (j in 1:ncol(emat)){
            e <- emat[,j]
            M[e[1], e[2]] <- M[e[2], e[1]] <- 1
        }
    }
    M
}

#' @export
#' @rdname coerce_models
as_emat2elist <- function(emat){
    if (ncol(emat) > 0)
        split(t.default(emat), 1:(ncol(emat)))
    else
        list()
}

#' @export
#' @rdname coerce_models
as_elist2emat <- function(elist){
    find_unique_rows <- function(ed){
        d <- max(ed)
        ii <- ed[,1] * 10^(nchar(d)+1) + ed[,2]
        ed[!duplicated(ii),]
    }
    
    idx <- sapply(elist, length) > 1
    elist <- elist[idx]
    
    if (length(elist) > 0){
        ed <- lapply(elist, combn_prim, 2)
        ed <- ed |> do.call(cbind, args=_)
        ed <- t.default(ed)
        b <- ed[, 1] > ed[, 2]
        ed[b,] <- ed[b, 2:1]
        ed <- find_unique_rows(ed)
        emat <- t.default(ed)
    } else
        emat <- matrix(NA, ncol=0, nrow=2)
    emat
}


#' @export
#' @rdname coerce_models
as_glist2emat <- function(glist){
  if (!inherits(glist, "list")) stop("'glist' must be a list\n")
  if (length(glist) == 0)
    return(matrix(NA, nrow=2, ncol=0))
  glist <- glist[sapply(glist, function(x) length(x)>1)]
  if (length(glist) == 0)
    return(matrix(NA, nrow=2, ncol=0))

  g2 <- lapply(glist, function(x) names2pairs(x))
  g2 <- unlist(g2, recursive=FALSE)
  g2 <- do.call(rbind, g2)
  t.default(g2[!duplicated(g2),])
}

#' @export
#' @rdname coerce_models
as_glist2cq <- function(glist){
    as_glist2emat(glist) |> as_emat2cq()
}


#' @export
#' @rdname coerce_models
as_glist2graph <- function(glist, d){
    as_emat2graph(as_glist2emat(glist), d=d)
}

#' @export
#' @rdname coerce_models
as_glist2igraph <- function(glist, d){
    g <- as_glist2graph(glist, d=d)
    as(g, "igraph")
}

#' @export
#' @rdname coerce_models
as_emat2graph  <- function(emat, d){
    as(as_emat2amat(emat, d=d), "igraph")
}

#' @export
#' @rdname coerce_models
as_emat2igraph <- function(emat, d){
    ## cat("as_emat2igraph\n")
    ## print(emat)
    igraph::make_undirected_graph(emat, d)
}

#' @export
#' @rdname coerce_models
as_amat2emat <- function(amat, eps=1e-4){
    amat[lower.tri(amat, diag=TRUE)] <- 0
    e <- which(abs(amat) > eps, arr.ind = T)
    rownames(e) <- NULL
    colnames(e) <- NULL
    t.default(e)
}

#' @export
#' @rdname coerce_models
as_emat2glist <- function(emat){
    if (ncol(emat)==0) return(list())
    gRbase::colmat2list(emat)
}

#' @export
#' @rdname coerce_models
as_glist2out_edges <- function(glist){
    u <- gRbase::ug(glist)
    n <- gRbase::nonEdgeList(u)
    do.call(rbind, lapply(n, as.numeric))
}

#' @export
#' @rdname coerce_models
as_K2amat <- function(K, eps=1e-4){
    MM <- zapsmall(K)
    MM <- 1 * (abs(MM) > eps)
    diag(MM) <- 0
    MM
}

#' @export
#' @rdname coerce_models
as_K2graph <- function(K){
    as(as_K2amat(zapsmall(K)), "igraph")
}

#' @export
#' @rdname coerce_models
as_sparse <- function(K){
    as(zapsmall(K), "dgCMatrix")}



## .form2emat <- function(form){
##     if (is.matrix(form) && nrow(form) == 2) return(form)  ## FIXME
##     else
##         if (length(form) > 0){
##             em <- do.call(rbind, lapply(form, function(g)
##                 if (length(g) > 1) t.default(combn_prim(g, 2))))
##             em <- t.default(unique(em))
##         } else
##             em  <- matrix(nrow=2, ncol=0)
##     em
## }

