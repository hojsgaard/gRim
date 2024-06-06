#' ---
#' title: Table 2 in paper
#' author: Søren Højsgaard and Steffen Lauritzen
#' date: Created "`r date()`"
#' output:
#'   html_document:
#'     toc: true
#' ---

##
## Tabel 2: covips1 og ncd1 (grids) ellers som nu.
## 

##+ TABLE 2

##+ Settings
##' # Settings

## rm(list = ls())
if (!exists("local_bench_settings"))
    local_bench_settings <- "bench_settings.r"
source(local_bench_settings)

tstart   <- Sys.time()
RDFILE_out <- paste0(.RES_DIR, "/", RDFILE2, ".RData")

##+ Create design for experiments
##' # Create design for experiments

grid_size <-
    matrix(c(
        12, 8,
        12, 16,
        24, 16,
        24, 32,
        48, 32,
        48, 64
    ), byrow=T, nc=2)

grid_size <-
    matrix(c(
        20, 25,
        40, 25,
        40, 50,
        80, 50
    ), byrow=T, nc=2)




grid_size <- round(grid_size / GRID_SCALING)
grid_size <- cbind(grid_size, apply(grid_size, 1, FUN=prod))
colnames(grid_size) <- c("nr", "nc", "nvar")

md_grid <- expand.grid(method = c("cov", "ncd")[1:2],
                       ver    = 1,
                       rep    = 1,
                       marg   = "edge",
                       dat    = c("n01", "pro"),
                       bench  = "B2",
                       stringsAsFactors = FALSE)

md_grid

design <- merge(md_grid, grid_size)

## stringr::str_pad(design$nvar, width="3", pad="0")
## sprintf("%.3f", v)

design_tmp <- design
design_tmp$nvar <- stringr::str_pad(design$nvar, width="4", pad="0")
nms <-
    design_tmp |>
    with(as.character(interaction(rep, method, marg, dat, nvar, bench, sep="_")))
rm(design_tmp)

row.names(design) <- nms
design

des_lst <- split_byrow(design)

##+ Create list of function calls, one for each row in the design
##' # Create list of function calls, one for each row in the design

arg_lst <- lapply(des_lst, function(r) {
    dat  <- data_lst[[r$dat, exact=FALSE]]
    da   <- dat[seq_len(r$nvar), seq_len(r$nvar)]
    g    <- as.numeric(r[c("nr", "nc")])
    gl   <- model_rectangular_grid(g, type="glist")
    list(S=da, formula=gl, method=r$meth, ver=r$ver)
})

fn_lst <- lapply(arg_lst, function(lsti){
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
## res_lst

##+ Create dataframe with result
##' # Create dataframe with result

raw <-
    res_lst |>
    lapply(summary) |>
    do.call(rbind, args=_)
raw$engine <- NULL
raw

raw$method <- NULL
result <- cbind(design, time=raw$time)  ## Requires things are in the right order

result <-
    result |>
    group_by(dat, marg, method, nr, nc, nvar) |>
    summarize(time=median(time)) |> print(n=30)

kable(result, format="pipe")


result2_wide <-
    result  |>
    arrange(dat, nvar, method)
result2_wide

result3_wide <- 
    result2_wide |> 
    pivot_wider(id_cols = c(nvar),
                names_from = c(dat, method),
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

