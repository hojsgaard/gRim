#' ---
#' title: Table 3 in paper
#' author: Søren Højsgaard and Steffen Lauritzen
#' date: Created "`r date()`"
#' output:
#'   html_document:
#'     toc: true
#' ---

##+ TABLE 3

##+ Settings
##' # Settings

## rm(list = ls())
## source(local_bench_settings)
if (!exists("local_bench_settings"))
    local_bench_settings <- "bench_settings.r"
source(local_bench_settings)
tstart   <- Sys.time()
RDFILE_out <- paste0(.RES_DIR, "/", RDFILE3, ".RData")


##+ Create design for experiments
##' # Create design for experiments

design<-
    expand.grid(
        method = c("cov", "ncd")[1:2],
        ## ver    = c(0),
        rep    = 1:TREE_NMOD,        
        marg   = c("edge"),
        dat    = c("pro", "n01"),
        prob   = TREE_PROB,
        nvar   = TREE_USE,
        eng    = ENG,
        bench  = "B3",
        stringsAsFactors = FALSE
    )

## design <- design[design$nvar==2000 & design$prob==0.01 & design$rep==5, ]

design <-
    design |>
    transform(Enparm = nvar + nvar * (nvar - 1) * prob / 2)
nms <-
    design |>
    with(as.character(interaction(rep, method, marg, dat, prob, nvar, bench, sep="_")))
## nms
nms <- gsub("_0.01_", "_0.010_", nms)

row.names(design) <- nms
design
design |> dim()

## Exclude too large configurations
design <-
    design |> filter(!(
        (prob == 0.005 & nvar %in% tail(TREE_USE,1)) |
        (prob == 0.010 & nvar %in% tail(TREE_USE,1)) 
    ))

## design |> xtabs(~nvar+prob, data=_)

## design_small <-
    ## design |> filter(nvar < tail(TREE_USE,1))

## design_large <-
    ## design |> filter(nvar == tail(TREE_USE,1))

## CHANGE THIS
## design <- design_small
## design <- design_large

design

design |> xtabs(~nvar+prob, data=_)

des_lst <- split_byrow(design)

all.equal(names(des_lst), rownames(design))

##+ Create list of function calls, one for each row in the design
##' # Create list of function calls, one for each row in the design

arg_lst <- lapply(des_lst, function(r){
    dat  <- data_lst[[r$dat, exact=FALSE]]
    da   <- dat[seq_len(r$n), seq_len(r$n)]
    set.seed(2022+r$rep)
    gl <- model_random_tree(r$nvar, prob=r$prob, type="glist")
    list(S=da, formula=gl, method=r$meth)
})

fn_lst <- lapply(arg_lst, function(lsti){
    set_default(fit_ggm, c(lsti, global_args))
})

##+ Fit models
##' # Fit models

t0 <- Sys.time()
res_lst <- 
    mcmapply(
        function(m, n){
            cat("stratum: ", n, "\n");
            out <- do.call(m, list())
            cat("stratum: ", n, "- done\n");            
            out
        }, fn_lst, names(fn_lst), SIMPLIFY = FALSE)
fit_time <- Sys.time() - t0

##+ Create dataframe with result
##' # Create dataframe with result

raw <-
    res_lst |>
    lapply(summary) |>
    do.call(rbind, args=_)
raw$engine <- NULL
raw

### Milisecs to secs
## raw$time <- raw$time / 1000

## raw$method <- NULL
result <- cbind(design, time=raw$time)  ## Requires things are in the right order

result <- result |>
    arrange(prob, marg, dat, nvar, method, time)
head(result, 100)

result <-
    result |>
    group_by(prob, marg, dat, nvar, method) |>
    summarize(time=median(time)) |> print(n=30)

kable(result, format="pipe")


result2_wide <-
    result  |>
    arrange(prob, nvar, dat, method)
result2_wide

result2_wide <- 
    result2_wide |> 
    pivot_wider(id_cols = c(marg, nvar, prob),
                names_from = c(method, dat),
                values_from = time)
result2_wide

kable(result2_wide, format="pipe")


#' save results
out <- list(sum=result, raw=raw, design=design)
save(out, file=RDFILE_out)

##+ Timing and saving result
##' # Timing and saving result

fit_time
script_time <- Sys.time() - tstart
script_time
sessionInfo()

