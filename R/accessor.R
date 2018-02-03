




getmi <- function(object, name=c("CGstats", "cgstats", "SSD", "ssd", "SS", "ss", "center")){
    
    switch(name,
           "CGstats"=,
           "cgstats"={object$datainfo$CGstats},
           "SSD"=,
           "ssd"={object$datainfo$CGstats$SSD},
           "SS"=,
           "ss"={object$datainfo$CGstats$SS},
           "center"               ={object$datainfo$CGstats$center},
           "disc.names"           =object$datainfo$disc.names,
           "cont.names"           =object$datainfo$cont.names,
           "glist"                =object$glist,
           "varNames"             =object$varNames,
           "S"                    =object$datainfo$S,
           "n"                    =object$datainfo$n,
           "data"                 =object$datainfo$data,
           "dimension"            =object$fitinfo$dimension,
           "dev"                  =object$fitinfo$dev,
           "logL"                 =object$fitinfo$logL,
           "aic"                  =object$fitinfo$aic,
           "bic"                  =object$fitinfo$bic,
           "isGraphical"          =object$properties["isg"],
           "isDecomposable"       =object$properties["issd"]
           )
}

.glist <- function(object){
    if (inherits(object, "iModel"))
        getmi(object, "glist")
    else if (inherits(object, "matrix"))
        getCliques(object)
    else stop("Do not know what to do\n")
}

.amat <- function(object, vn = NULL, result = "matrix") {
    if (inherits(object, c("list", "formula"))){
        glist <- rhsf2list(object)
        if (is.null(vn))
            vn <- unique.default(c(glist, recursive=TRUE))
        .glist2amat(glist, vn, result)
    } else if (inherits(object, "iModel")) {
        .glist2amat(.glist(object), vn, result)        
    }
    else stop("Do not know what to do\n")       
}


.glist2amat <- function (glist, vn = NULL, result = "matrix") {
    #cat("caller of .glist2amat: ", deparse(sys.calls()[[sys.nframe()-1]]), "\n")
    glist2adjMAT(glist, vn = vn, result = result)
}

.as_amat <- function(x, vn=NULL){
    ##cat("caller of .as_amat: ", deparse(sys.calls()[[sys.nframe()-1]]), "\n")
    if (inherits(x, "list")){
        if (is.null(vn)) vn <- unique.default(c(x, recursive=TRUE))
        .glist2amat(x, vn=vn)
    } else if (inherits(x, "iModel"))
        .glist2amat(getmi(x, "glist"))      
}
