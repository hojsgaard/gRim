#' ---
#' title: Table 1 in paper
#' author: Søren Højsgaard and Steffen Lauritzen
#' date: Created "`r date()`"
#' output:
#'   html_document:
#'     toc: true
#' ---

##
## Tabel 1: covips1 og conips1 og ncd1 for 100 var og cov-con med edge-clique
## 

##+ TABLE 1
## rm(list = ls())
if (!exists("local_bench_settings"))
    local_bench_settings <- "bench_settings.r"
source(local_bench_settings)

tstart   <- Sys.time()
RDFILE_out <- paste0(.RES_DIR, "/", RDFILE1, ".RData")

#' ## Settings
#+ settings
NVAR_VEC
PROB_VEC
NVAR_FIX
NMOD
EPS

##+ Create design for experiments
##' # Create design for experiments

design <-
    expand.grid(
        method = c("cov", "con", "ncd")[1:3],
        ## ver    = c(1),
        marg   = c("cliq", "edge"),
        rep    = 1:NMOD,
        dat    = c("n01", "pro"),
        prob   = PROB_VEC,
        nvar   = NVAR_FIX,
        bench  = "B1",
        stringsAsFactors = FALSE
    )
names(design)

design <-
    design |>
    transform(Enparm = nvar + nvar * (nvar - 1) * prob / 2)
nms <-
    design |>
    with(as.character(interaction(rep, method, marg, dat, prob, nvar, bench, sep="_")))
row.names(design) <- nms
design |> head(30)

design <-
    design |> filter(!(method=="ncd" & marg=="cliq"))
design


##+ Create list of function calls, one for each row in the design
##' # Create list of function calls, one for each row in the design

des_lst <- split_byrow(design)
arg_lst <- lapply(des_lst, function(r) {
    dat  <- data_lst[[r$dat, exact=FALSE]]
    da   <- dat[seq_len(r$n), seq_len(r$n)]
    set.seed(2022 + r$rep)
    em   <- emat_random_model(1:r$nvar, prob = r$prob)
    if (identical(r$marg, "cliq"))
        gl <- as_emat2cq(em)
    else
        gl   <- as_emat2glist(em)
    list(S=da, formula=gl, method=r$meth)
})

fn_lst <- lapply(arg_lst, function(lsti) {
    set_default(fit_ggm, c(lsti, global_args))
})

##+ Fit models
##' # Fit models

t0 <- Sys.time()
res_lst <- 
    mcmapply(
        function(m, n) {
            cat("stratum: ", n, "\n");
            out <- do.call(m, list())
            cat("stratum: ", n, "- done\n");
            out
        }, fn_lst, names(fn_lst), SIMPLIFY = FALSE)
fit_time <- Sys.time() - t0
res_lst


## FIXME Save fn_lst and res_lst ???

##+ Create dataframe with result
##' # Create dataframe with result

raw <-
    res_lst |>
    lapply(summary) |>
    do.call(rbind, args=_)
raw$eng <- NULL
raw

raw$method <- NULL
result <- cbind(design, time=raw$time)  ## Requires things are in the right order

result <- arrange(result, dat, method)
result

result <-
    result |>
    group_by(dat, prob, marg, method, nvar) |>
    summarize(time=median(time)) |> print(n=30)

kable(result, format="pipe")


result2_wide <-
    result  |>
    arrange(dat, marg, method)
result2_wide

result3_wide <- 
    result2_wide |> 
    pivot_wider(id_cols = c(prob),
                names_from = c(method, marg, dat),
                values_from = time)
result3_wide

kable(result3_wide, format="pipe")



#' save results

out <- list(sum=result, raw=raw, design=design)
save(out, file=RDFILE_out)

##+ Timing and saving result
##' # Timing and saving result

fit_time
script_time <- Sys.time() - tstart
script_time
sessionInfo()


