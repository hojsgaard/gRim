###
###
### NEW - october 2015
###
###

##.getmi <- function(object, name=c("CGstats", "cgstats", "SSD", "ssd", "SS", "ss", "center")){
##    
##    switch(name,
##           "CGstats"=, "cgstats"={object$datainfo$CGstats},
##           "SSD"=,"ssd"={object$datainfo$CGstats$SSD},
##           "SS"=,"ss"={object$datainfo$CGstats$SS},
##           "center"={object$datainfo$CGstats$center},
##           )
##}
##

.fit222 <- function(object, method="general", details=0, eps.parm=1e-6, maxit=100,...){
  method <- match.arg(method, c("general", "stephalving"))

  ## Saturated model
  n.obs   <- getmi(object, "CGstats")$n.obs
  NN      <- getmi(object, "CGstats")$N
  Sig.sat <- getmi(object, "CGstats")$SSD/NN
  qq        <- length(getmi(object, "CGstats")$cont.names)
  logL.sat  <- sum(n.obs*log(n.obs/NN)) - NN*qq/2 * log(2*pi) - NN/2 * log(det(Sig.sat)) - NN*qq/2

  
  ## Independence model
  i.model    <- loglin(n.obs, as.list(seq_along(dim(n.obs))),iter=1, print=FALSE, fit=TRUE)
  grand.mean <- rowSumsPrim(colwiseProd(n.obs/NN, getmi(object, "center")))
  Sig.ind    <- (getmi(object, "SS") - NN*grand.mean %*% t(grand.mean))/NN
  logL.ind   <- sum(n.obs*log(i.model$fit/NN))- NN*qq/2 * log(2*pi) - NN/2 * sum(log(diag(Sig.ind))) - NN*qq/2

  ans           <- .fitmModel222(object, method=method, details=details, eps.parm=eps.parm, maxit=maxit, ...)
  ans$dimension <- .mmod_dimension(object$modelinfo$dlq, object$datainfo)
  
  ans$ideviance <-  ans$logL - logL.ind
  ans$dev       <-  logL.sat - ans$logL
  ans$aic       <-  -2 * ans$logL + 2*ans$dimension['mod.dim']
  ans$bic       <-  -2 * ans$logL + log( NN )*ans$dimension['mod.dim']

  object$fitinfo  <- ans
  object$isFitted <- TRUE
  return(object)
}

.weakMarginalData <- function(object, type="ghk", details=0){

  .infoPrint(details,1,"fit.mModel: Finding weak (empirical) marginals for each generator\n")

  dgen.list <- object$modelinfo$glist.num.disc
  cgen.list <- object$modelinfo$glist.num.cont
  CGstats   <- getmi(object, "CGstats")
  
  ans <- vector("list", length(dgen.list))
  for ( i  in seq_along(dgen.list)){
    EE.mm <- weakMarginalData(CGstats, disc=dgen.list[[ i ]], cont=cgen.list[[ i ]],
                              type=type, details=details)
    ans[[ i ]] <- EE.mm
  }
  ans
}


.fitmModel222 <- function(object, method="general", details=0, eps.parm=1e-6, eps.logL=1e-6, maxit=100,...){
  t.start <- proc.time()
  ##cgs<<-getmi(object, "CGstats")
  #cat("+++ CGstats:\n"); print(names(getmi(object, "CGstats")))
  ### FIXME: Here SSD and SS are available
 
### Generator lists (numerical values)
  dgen.list <- object$modelinfo$glist.num.disc
  cgen.list <- object$modelinfo$glist.num.cont
  
### Find weak marginal model for each generator
  ghk.obs.list <- .weakMarginalData(object, type="ghk", details=details)
  pms.obs.list <- lapply(ghk.obs.list, ghk2pmsParms_)
  
  ### Generate initial parameter value
  .infoPrint(details,1, "fit.mModel: Creating initial values\n")
  Mparms <- .CGstats2initpms( getmi(object, "CGstats") )[c(1:4,8)]
  Cparms <- pms2ghkParms(Mparms)

  #' init <<- list(Cparms=Cparms, Mparms=Mparms,
  #'               dgen.list=dgen.list, cgen.list=cgen.list,
  #'               ghk.obs.list=ghk.obs.list, pms.obs.list=pms.obs.list)
  
  
### Iterate to maximize logL
  .infoPrint(details,1, "fit.mModel: Calling .mModel_iterate\n")
  ## .mModel_iterate returns an update of Cparms + a bit more
  ans <- .mModel_iterate222(list(Mparms=Mparms, Cparms=Cparms),
                            dgen.list, cgen.list,
                            ghk.obs.list, pms.obs.list,
                            getmi(object, "CGstats")[c("n.obs","center","SSD")], 
                            method, eps.parm, eps.logL, maxit, details)

  ans$ghk.obs.list <- ghk.obs.list
  ans$init.parms   <- Mparms
  class(ans)       <- "MIfit"
  
  t.use <- proc.time()-t.start #; cat("Fitting time:",t.use[1],"\n")
  ans
}

.mModel_iterate222 <- function(MCparms, dgen.list, cgen.list,
                               ghk.obs.list, pms.obs.list, CGstats,
                               method="general",
                               eps.parm=1.0e-4, eps.logL=1e-6, maxit=100, details=1){
    
    prev.Mparms <- MCparms$Mparms
    init.logL   <- prev.logL <- .mModel_logLpms(CGstats, MCparms$Mparms)

    .infoPrint(details,3, cat(sprintf("Initial logL %f\n", init.logL)))

    itcount    <- 1
    scale      <- 1
    
    if (maxit>0){
        repeat {
            zzz       <- .outerloop222(MCparms, dgen.list, cgen.list,
                                       ghk.obs.list, pms.obs.list, CGstats, scale,
                                       method, prev.logL, itcount, details)
            curr.logL <- zzz$logL
            d.logL    <- zzz$d.logL
            d.parms   <- zzz$d.parms
            
            .infoPrint(details,1,
                       cat(sprintf(".mModel_iterate maxit=%i, curr.logL=%14.6f, d.logL=%12.6f, d.parms=%12.8f\n",
                                   maxit, curr.logL, d.logL, d.parms)))
            
            it.exceed.crit <- itcount >=maxit
            neg.d.logL     <- d.logL < -1e-6

            #print(d.logL, digits=20)
            if (!neg.d.logL & d.logL < eps.logL){
                break() # We are done
            } else {
                if (neg.d.logL){
                    cat(sprintf("Fitting method=%s; logL failed to increase (d.logL=%f) - fit may be questionable\n",
                                method,d.logL))
                    break()
                } else {
                    if (it.exceed.crit){
                        cat(sprintf("Fitting method=%s; Maximum number of iterations=%i exceeded - fit may be questionable\n",
                                    method,itcount))
                        break()
                    }
                }
            }
            
            #' Mparms      <- zzz$Mparms
            #' Cparms      <- zzz$Cparms
            MCparms     <- list(Mparms=zzz$Mparms, Cparms=zzz$Cparms)
            prev.Mparms <- zzz$Mparms
            prev.logL   <- curr.logL
            itcount     <- itcount + 1
        }
    }
    list(parms=zzz$Cparms, logL=curr.logL, init.logL=init.logL)
}


.outerloop222 <- function(MCparms, dgen.list, cgen.list,
                          ghk.obs.list, pms.obs.list, CGstats, scale,
                          method, logL, itcount, details){


    .infoPrint(details,4, "calling outerloop\n")
    
    prev.Mparms <- MCparms$Mparms
    prev.logL   <- logL
    logL.fail   <- 0

    ##pmp.new <<- prev.Mparms
    if (method=="general")
        .innerloop <- .standard.innerloop222
    else
        .innerloop <- .stephalving.innerloop222

    ## Do one iteration over all generators
    i <- 1
    while( i <= length( dgen.list )){
        dgen.idx    <- dgen.list[[ i ]]
        cgen.idx    <- cgen.list[[ i ]]
        ghk.obs     <- ghk.obs.list[[ i ]]
        pms.obs     <- pms.obs.list[[ i ]]
        pms.fit     <- weakMarginalModel(MCparms$Mparms, disc=dgen.idx, cont=cgen.idx, type="pms", details=details)
        ghk.fit     <- pms2ghkParms(pms.fit)
        upd         <- .innerloop(MCparms, dgen.idx, cgen.idx,
                                  ghk.obs, pms.obs, ghk.fit, pms.fit, CGstats, scale, prev.logL, details)
        #' Mparms    <- upd$Mparms
        #' Cparms    <- upd$Cparms
        MCparms   <- list(Mparms=upd$Mparms, Cparms=upd$Cparms)
        #' curr.logL <- upd$curr.logL              # FIXME What is this
        #' logL.fail <- logL.fail + upd$logL.fail  # FIXME What is this
        i <- i + 1
    }

###    mp.new <<- upd$Mparms
   ## stop()
    
    curr.logL <- .mModel_logLpms( CGstats, upd$Mparms )
    d.logL    <- curr.logL - prev.logL
    d.parms   <- .mModel_parmdiff222(upd$Mparms, prev.Mparms)
    
    .infoPrint(details,3,
               cat(sprintf("outerloop (%2d): logL %16.10f, d.logL: %16.10f d.parms: %16.10f logL.fail: %f\n",
                           itcount, curr.logL, d.logL, d.parms, logL.fail)))
    
    list(Mparms=upd$Mparms, Cparms=upd$Cparms, logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail, d.parms=d.parms)
}



.standard.innerloop222 <- function(MCparms, dgen.idx, cgen.idx,
                                   ghk.obs, pms.obs, ghk.fit, pms.fit, CGstats, scale, prev.logL, details){
    
    new.Cparms <- .update.ghkParms222(MCparms$Cparms, dgen.idx, cgen.idx,
                                      ghk.obs, pms.obs, ghk.fit, pms.fit, scale, details=details)

    ##ncp.new <<- new.Cparms
    ##new.Mparms <- ghk2pmsParms(new.Cparms)
    new.Mparms <- ghk2pmsParms_(new.Cparms) 
    
    curr.logL  <- d.logL <- d.parms <- NA
    logL.fail  <- as.numeric(d.logL < 0)
    
    .infoPrint(details,4,
               cat(sprintf(".std.innerloop(%4.2f): AA=%10s,  curr.logL=%16.10f -2logL=%16.10f d.logL=%16.10f d.parms=%8.6f \n",
                           scale, .toString(c("{",dgen.idx,"|", cgen.idx,"}")),
                           curr.logL, -2*curr.logL, d.logL,
                           d.parms)))
    
    list(Mparms=new.Mparms, Cparms=new.Cparms, curr.logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail,
         maxinner.code=NA, step.code=NA)
}


.stephalving.innerloop222 <- function(MCparms, dgen.idx, cgen.idx,
                                      ghk.obs, pms.obs, ghk.fit, pms.fit, CGstats, scale, prev.logL, details){
    
    .infoPrint(details,10, "innerloop: finding (model) weak marginals for a generator\n")
    
    prev.Mparms   <- MCparms$Mparms
    innercount    <-  1
    maxinner      <-  5
    neg.eps       <-  -1e-4
    good.Mparms   <- MCparms$Mparms
    good.Cparms   <- MCparms$Cparms
    step.code     <- 0
    maxinner.code <- 0
    
    d.logL        <- -99999
    curr.logL     <- prev.logL
    
    repeat{
        new.Cparms <- .update.ghkParms222(good.Cparms, dgen.idx, cgen.idx,
                                          ghk.obs, pms.obs, ghk.fit, pms.fit, scale, details=details)
        ##new.Mparms <- ghk2pmsParms(new.Cparms)
        new.Mparms <- ghk2pmsParms_(new.Cparms)

        ## FIXME: 2016-11-06: Line below
        ##curr.logL  <- .mModel_logLpms222( CGstats, new.Mparms )
        ## FIXME: 2016-11-06: changed to (because then fun does not exist):
        curr.logL  <- .mModel_logLpms( CGstats, new.Mparms )
        d.logL     <- curr.logL - prev.logL
        d.parms    <- .mModel_parmdiff222( new.Mparms, prev.Mparms )
        min.eigen  <- min(eigen(new.Mparms$Sigma)$values)
        
        .infoPrint(details,1,
                   cat(sprintf(".steph.innerloop(%4.2f): AA=%10s,  curr.logL=%16.10f -2logL=%16.10f d.logL=%16.10f d.parms=%8.6f \n",
                               scale, .toString(c("{",dgen.idx,"|", cgen.idx,"}")), curr.logL, -2*curr.logL, d.logL,
                               d.parms)))
        
        if ((d.logL<neg.eps | min.eigen<0) & innercount < maxinner){
            scale       <- scale / 2
            innercount  <- innercount + 1
            step.code   <- 1
        } else {
            if (innercount==maxinner){
                Cparms    <- good.Cparms
                Mparms    <- good.Mparms
                curr.logL <- .mModel_logLpms(CGstats, Mparms)
                maxinner.code <- 1
                .infoPrint(details, 4, cat(sprintf("stephalving failed; restoring original parameters; logL: %10.4f\n", curr.logL)))
            } else {
                Cparms  <- new.Cparms
                Mparms  <- new.Mparms
            }
            break
        }
    }
    
    logL.fail <- as.numeric(d.logL < 0)
    list(Mparms=Mparms, Cparms=Cparms, curr.logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail,
                maxinner.code=maxinner.code, step.code=step.code)
}

## CGstats argument never used.
.update.ghkParms222 <- function(Cparms, dgen.idx, cgen.idx, ghk.obs, pms.obs, ghk.fit, pms.fit, scale, details=0) {

    #' cat(".update.ghkParms: Cparms:\n"); print(Cparms)
    ## cat("ghk.fit:\n"); print(ghk.fit)
    
    g.idx <- 1L
    h.idx <- 2L
    K.idx <- 3L

    #normalize.ghkParms <- .normalize.ghkParms #normalize_ghkParms_
    normalize.ghkParms <- normalize_ghkParms_
#    ghk2pmsParmsFUN    <- ghk2pmsParms_
    
    .infoPrint(details,5, cat(sprintf(".update.ghkParms: A=%8s\n",
                                      .toString(c("{",dgen.idx,"|", cgen.idx,"}")))))
    #' gt <- .genType(dgen.idx, cgen.idx) ## ; cat("generator", gt, "\n")
    #' print(gt)
    #' print(ghk.fit$gentype)
    gt <- ghk.fit$gentype
    ##gt <- Cparms$gentype

    ##print(gt)
    d.parms.crit <- 0.00001
    
    if (details>=5){
        cat("PRE UPDATED marginal OBSERVED // FITTED values - moment form\n")
        print(rbind(.as.matrix(ghk2pmsParms_( ghk.obs )),.as.matrix(ghk2pmsParms_(ghk.fit))))
    }
    
    ##  cat(".update.ghkParms - calling .mModel_parmdiff\n")
    marg.d.parms <- .mModel_parmdiff222(pms.fit, pms.obs)
    .infoPrint(details,5, cat(sprintf("PARMDIF=%f\n", marg.d.parms)))
    
    if (marg.d.parms <= d.parms.crit){
        .infoPrint(details, 5, cat(sprintf("Not updating generator\n")))
        res <- Cparms
    } else {
        if (details>=5){
            cat("PRE UPDATE Mparms:\n");
            print(.MIparms2matrix(ghk2pmsParms_( Cparms )))
            cat("PRE UPDATED marginal OBSERVED // FITTED values - canonical form\n")
            print(rbind(.as.matrix((ghk.obs)),.as.matrix((ghk.fit))))
            cat("PRE UPDATE Cparms:\n");
            print(.MIparms2matrix((Cparms)))
        }
        
        switch(gt,
               "discrete"={
                   ##cat("Cparms//Mparms BEFORE update:\n");  print(Cparms)
                   upd.g    <- scale*(ghk.obs[[ g.idx ]] - ghk.fit[[ g.idx ]])
                   g.new    <- tableOp2(Cparms[[ g.idx ]], upd.g, `+`, restore=TRUE)
                   res <- list(g=g.new, h=Cparms[[ h.idx ]], K=Cparms[[ K.idx ]], gentype="discrete")
                   ##r3 <<- res
                   ##cat("Cparms//Mparms AFTER update:\n");  print(res)
               },
               "continuous"={
                   h.new   <- Cparms[[ h.idx ]]
                   upd.h   <- scale * (ghk.obs[[ h.idx ]]-ghk.fit[[ h.idx ]])
                   for (jj in 1:ncol(h.new))
                       h.new[cgen.idx, jj] <- Cparms[[ h.idx ]][cgen.idx, jj, drop=FALSE] + upd.h
                   
                   upd.k   <- scale * (ghk.obs[[ K.idx ]] - ghk.fit[[ K.idx ]])
                   K.new   <- Cparms[[ K.idx ]]
                   K.new[cgen.idx, cgen.idx] <- K.new[cgen.idx, cgen.idx] + upd.k
                   ## cat("cont: upd.h:\n"); print(cbind(ghk.obs[[ h.idx ]], ghk.fit[[ h.idx ]], upd.h))
                   ## cat("cont: upd.k:\n"); print(cbind(ghk.obs[[ K.idx ]], ghk.fit[[ K.idx ]], upd.k))
                   res <- list(g=Cparms[[ g.idx ]], h=h.new, K=K.new, gentype="continuous")
               },
               "mixed"={
                   ## g update:
                   upd.g    <- scale * (ghk.obs[[g.idx]] - ghk.fit[[g.idx]])
                   g.new    <- tableOp2(Cparms[[g.idx]], upd.g, `+`, restore=TRUE)                 

                   ##cat("upd.g:\n"); print(t(round(cbind(ghk.obs[["g"]], ghk.fit[["g"]], upd.g),4)))
                   ## K update:
                   upd.k   <- scale * (ghk.obs[[K.idx]] - ghk.fit[[K.idx]])
                   K.new   <- Cparms[[K.idx]]
                   K.new[cgen.idx,cgen.idx] <- K.new[cgen.idx,cgen.idx] + upd.k
                   ##cat("upd.k:\n"); print(round(cbind(ghk.obs[[ K.idx ]], ghk.fit[[ K.idx ]], upd.k),4))
                   h.new   <- Cparms[[h.idx]]
                   upd.h   <- scale * (ghk.obs[[h.idx]]-ghk.fit[[h.idx]])
                   
                   em       <- ghk.fit[['jia.mat']]
                   Cparms.h <- Cparms[[ h.idx ]]
                   for (jj in 1:ncol(em))
                       h.new[cgen.idx,em[,jj]] <- Cparms.h[cgen.idx,em[,jj],drop=FALSE] + upd.h[,jj]
                   ##cat("upd.h:\n"); print(round(cbind(ghk.obs[[ h.idx ]], ghk.fit[[ h.idx ]], upd.h),4))
                   ##max.chg <- c(max(abs(upd.g)),max(abs(upd.h)),max(abs(upd.k)))
                   res <- list(g=g.new, h=h.new, K=K.new, gentype="mixed")
               })
        ##res <<- res
        res <- normalize_ghkParms_(res)
        #res <- .normalize.ghkParms(res)
        res <- c(res[1:3], Cparms[-(1:3)])      
    }
    
    if (details>=5){
        cat("POST UPDATE Cparms // Mparms:\n");
        MM   <- ghk2pmsParms_( Cparms )
        MM$p <- MM$p*MM$N
        RR   <- ghk2pmsParms_( res )
        RR$p <- RR$p*RR$N
        print(rbind(.as.matrix(res), .as.matrix(RR)))
    }
                                        #res$max.chg <- max.chg
    #' cparms.new <<- Cparms
    #' res.new <<- res
    res
}


.mModel_logLpms2 <- function(CGstats, Mparms){

  Sigma.inv <- solveSPD(Mparms[['Sigma']])
  n.i    <- as.numeric(CGstats[['n.obs']])
  N      <- sum(n.i)
  Q      <- nrow(CGstats[['center']])
  xxx    <- sum(n.i * log(Mparms[['p']])) - N * (Q * log(2*pi) + .logdet(Mparms[['Sigma']])) / 2
  x4     <- -  sum(CGstats[['SSD']] * Sigma.inv) / 2
  mu.dif <- CGstats[['center']] - Mparms[['mu']]
  quad   <- .vMMt(n.i, mu.dif)
  x5     <- - sum(Sigma.inv * quad) / 2
  xxx + x4 + x5
}


.mModel_parmdiff222 <- function(curr.Mparms, prev.Mparms){

    #' cat("curr.Mparms:---------------\n "); print(curr.Mparms)
    #' cat("prev.Mparms:---------------\n "); print(prev.Mparms)

    if (curr.Mparms[['gentype']]=="discrete"){
        N   <- prev.Mparms[['N']]
        cp  <- curr.Mparms[['p']]
        ppp <- as.numeric(N * abs((cp - prev.Mparms[['p']])) /sqrt((N * cp + 1)))
        ans <- max(ppp)
    } else {
        N   <- prev.Mparms[['N']]
        cp  <- curr.Mparms[['p']]
        #sss <- curr.Mparms[['Sigma']]
        #nr  <- nrow(sss)
        nr  <- nrow( curr.Mparms[['Sigma']] )
        iii <- 1+(nr+1)*((1:nr)-1)
        #ddd <- sss[iii] ## faster than diag(sss)
        ddd <- curr.Mparms[['Sigma']][iii] ## faster than diag(sss)
        ##     ppp <- as.numeric(N * abs((cp - prev.Mparms[[p.idx]])) / sqrt((N * cp + 1)))
        ##     mmm <- as.numeric(abs(curr.Mparms[[mu.idx]] - prev.Mparms[[mu.idx]])/sqrt(ddd))
        
        #' vv  <<- list(N=N, cp=cp, pp=prev.Mparms[['p']])
        #' print(vv)
        ppp <- c(N * abs((cp - prev.Mparms[['p']])) / sqrt((N * cp + 1)))
        mmm <- c(abs(curr.Mparms[['mu']] - prev.Mparms[['mu']])/sqrt(ddd))
        
        #' xxx <- abs(sss - prev.Mparms[['Sigma']])
        #' uuu <- xxx/sqrt(tcrossprod(ddd) + sss^2)

        xxx <- abs(curr.Mparms[['Sigma']] - prev.Mparms[['Sigma']])
        uuu <- xxx/sqrt(tcrossprod(ddd) + curr.Mparms[['Sigma']]^2)

        ans <- max(c(ppp,mmm,c(uuu)))
        ##  cat("max parm diff:", ans, "\n")
    }
    ans
}




#' .stephalving.innerloop222 <- function(Mparms, Cparms, dgen.idx, cgen.idx,
#'                                       ghk.obs, pms.obs, ghk.fit, pms.fit, CGstats, scale, prev.logL, details){
    
#'     .infoPrint(details,10, "innerloop: finding (model) weak marginals for a generator\n")
    
#'     prev.Mparms   <- Mparms
#'     innercount    <-  1
#'     maxinner      <-  5
#'     neg.eps       <-  -1e-4
#'     good.Mparms   <- Mparms
#'     good.Cparms   <- Cparms
#'     step.code     <- 0
#'     maxinner.code <- 0
    
#'     d.logL      <- -99999
#'     curr.logL   <- prev.logL
    
#'     repeat{
#'         new.Cparms <- .update.ghkParms222(good.Cparms, dgen.idx, cgen.idx,
#'                                           ghk.obs, pms.obs, ghk.fit, pms.fit, scale, details=details)
#'         ##new.Mparms <- ghk2pmsParms(new.Cparms)
#'         new.Mparms <- ghk2pmsParms_(new.Cparms)
        
#'         curr.logL  <- .mModel_logLpms222( CGstats, new.Mparms )
#'         d.logL     <- curr.logL - prev.logL
#'         d.parms    <- .mModel_parmdiff222( new.Mparms, prev.Mparms )
#'         min.eigen  <- min(eigen(new.Mparms$Sigma)$values)
        
#'         .infoPrint(details,1,
#'                    cat(sprintf(".steph.innerloop(%4.2f): AA=%10s,  curr.logL=%16.10f -2logL=%16.10f d.logL=%16.10f d.parms=%8.6f \n",
#'                                scale, .toString(c("{",dgen.idx,"|", cgen.idx,"}")), curr.logL, -2*curr.logL, d.logL,
#'                                d.parms)))
        
#'         if ((d.logL<neg.eps | min.eigen<0) & innercount < maxinner){
#'             scale       <- scale / 2
#'             innercount  <- innercount + 1
#'             step.code   <- 1
#'         } else {
#'             if (innercount==maxinner){
#'                 Cparms    <- good.Cparms
#'                 Mparms    <- good.Mparms
#'                 curr.logL <- .mModel_logLpms(CGstats, Mparms)
#'                 maxinner.code <- 1
#'                 .infoPrint(details, 4, cat(sprintf("stephalving failed; restoring original parameters; logL: %10.4f\n", curr.logL)))
#'             } else {
#'                 Cparms  <- new.Cparms
#'                 Mparms  <- new.Mparms
#'             }
#'             break
#'         }
#'     }
    
#'     logL.fail <- as.numeric(d.logL < 0)
#'     list(Mparms=Mparms, Cparms=Cparms, curr.logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail,
#'                 maxinner.code=maxinner.code, step.code=step.code)
#' }


#' .standard.innerloop222 <- function(Mparms, Cparms, dgen.idx, cgen.idx,
#'                                    ghk.obs, pms.obs, ghk.fit, pms.fit, CGstats, scale, prev.logL, details){
    
#'     new.Cparms <- .update.ghkParms222(Cparms, dgen.idx, cgen.idx,
#'                                       ghk.obs, pms.obs, ghk.fit, pms.fit, scale, details=details)
    
#'     ##new.Mparms <- ghk2pmsParms(new.Cparms)
#'     new.Mparms <- ghk2pmsParms_(new.Cparms) 
    
#'     curr.logL  <- d.logL <- d.parms <- NA
#'     logL.fail  <- as.numeric(d.logL < 0)
    
#'     .infoPrint(details,4,
#'                cat(sprintf(".std.innerloop(%4.2f): AA=%10s,  curr.logL=%16.10f -2logL=%16.10f d.logL=%16.10f d.parms=%8.6f \n",
#'                            scale, .toString(c("{",dgen.idx,"|", cgen.idx,"}")), curr.logL, -2*curr.logL, d.logL,
#'                            d.parms)))
    
#'     list(Mparms=new.Mparms, Cparms=new.Cparms, curr.logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail,
#'          maxinner.code=NA, step.code=NA)
#' }


#' .mModel_iterate222 <- function(Mparms, Cparms, dgen.list, cgen.list,
#'                                ghk.obs.list, pms.obs.list, CGstats,
#'                                method="general",
#'                                eps.parm=1.0e-4, eps.logL=1e-6, maxit=100, details=1){
    
#'     prev.Mparms <- Mparms
#'     init.logL   <- prev.logL <- .mModel_logLpms(CGstats, Mparms)

#'     .infoPrint(details,3, cat(sprintf("Initial logL %f\n", init.logL)))

#'     itcount    <- 1
#'     scale      <- 1
    
#'     if (maxit>0){
#'         repeat {
#'             zzz       <- .outerloop222(Mparms, Cparms, dgen.list, cgen.list,
#'                                        ghk.obs.list, pms.obs.list, CGstats, scale,
#'                                        method, prev.logL, itcount, details)
#'             curr.logL <- zzz$logL
#'             d.logL    <- zzz$d.logL
#'             d.parms   <- zzz$d.parms
            
#'             .infoPrint(details,1,
#'                        cat(sprintf(".mModel_iterate maxit=%i, curr.logL=%14.6f, d.logL=%12.6f, d.parms=%12.8f\n",
#'                                    maxit, curr.logL, d.logL, d.parms)))
            
#'             it.exceed.crit <- itcount >=maxit
#'             neg.d.logL     <- d.logL < 0
            
#'             if (!neg.d.logL & d.logL < eps.logL){
#'                 break() # We are done
#'             } else {
#'                 if (neg.d.logL){
#'                     cat(sprintf("Fitting method=%s; logL failed to increase (d.logL=%f) - fit may be questionable\n",
#'                                 method,d.logL))
#'                     break()
#'                 } else {
#'                     if (it.exceed.crit){
#'                         cat(sprintf("Fitting method=%s; Maximum number of iterations=%i exceeded - fit may be questionable\n",
#'                                     method,itcount))
#'                         break()
#'                     }
#'                 }
#'             }
            
#'             Mparms      <- zzz$Mparms
#'             Cparms      <- zzz$Cparms
#'             prev.Mparms <- Mparms
#'             prev.logL   <- curr.logL
#'             itcount     <- itcount + 1
#'         }
#'     }
#'     list(parms=Cparms, logL=curr.logL, init.logL=init.logL)
#' }


#' .outerloop222 <- function(Mparms, Cparms, dgen.list, cgen.list,
#'                           ghk.obs.list, pms.obs.list, CGstats, scale,
#'                           method, logL, itcount, details){

#'     .infoPrint(details,4, "calling outerloop\n")
    
#'     prev.Mparms <- Mparms
#'     prev.logL   <- logL
#'     logL.fail   <- 0
    
#'     if (method=="general")
#'         .innerloop <- .standard.innerloop222
#'     else
#'         .innerloop <- .stephalving.innerloop222
    
#'     i <- 1
#'     while( i <= length( dgen.list )){
#'         ##for ( i  in seq_along( dgen.list )){
#'         dgen.idx    <- dgen.list[[ i ]]
#'         cgen.idx    <- cgen.list[[ i ]]
#'         ghk.obs     <- ghk.obs.list[[ i ]]
#'         pms.obs     <- pms.obs.list[[ i ]]
#'         pms.fit     <- weakMarginalModel(Mparms, disc=dgen.idx, cont=cgen.idx, type="pms", details=details)
#'         ghk.fit     <- pms2ghkParms(pms.fit)
#'         zzz         <- .innerloop(Mparms, Cparms, dgen.idx, cgen.idx,
#'                                   ghk.obs, pms.obs, ghk.fit, pms.fit, CGstats, scale, prev.logL, details)
#'         Mparms    <- zzz$Mparms
#'         Cparms    <- zzz$Cparms
#'         curr.logL <- zzz$curr.logL
#'         logL.fail <- logL.fail + zzz$logL.fail
#'         i <- i + 1
#'     }
    
#'     curr.logL <- .mModel_logLpms(CGstats, Mparms)
#'     d.logL    <- curr.logL - prev.logL
#'     d.parms   <- .mModel_parmdiff222(Mparms, prev.Mparms)
    
#'     .infoPrint(details,3,
#'                cat(sprintf("outerloop (%2d): logL %16.10f, d.logL: %16.10f d.parms: %16.10f logL.fail: %f\n",
#'                            itcount, curr.logL, d.logL, d.parms, logL.fail)))
    
#'     list(Mparms=Mparms, Cparms=Cparms, logL=curr.logL, d.logL=d.logL, logL.fail=logL.fail, d.parms=d.parms)
#' }
