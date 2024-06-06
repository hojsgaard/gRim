##' -----------------------------------------------------------------
##' # Benchmarking for paper "On some algorithms for estimation in
##' # Gaussian graphical models""
##' 
##' # Copyright: Søren Højsgaard and Steffen Lauritzen
##'
##' Notice: Running these scripts on a 64 cores linux platform takes
##' about 34 hours. 
##' -----------------------------------------------------------------

##' design: (either of 1,2,3; i.e. 1 = the design in the paper, 2 = a small
##' (and fast) design and 3 = a larger (and slower) design)
design <- 1

##' table: Which table(s) to genrate ## (a subset of c(1, 2, 3))
table  <- 2

##' Global settings
EPS = 1e-3  


## #####################################################
## NEED NOT TOUCH ANYTHING BELOW HERE 
## #####################################################

##' ## SET VARIABLES TO BE USED IN SCRIPTS
use <- c("design_paper",   ## Settings generating data in paper
         "design_small",   ## Small settings for testing purposes
         "design_medium"   ## Larger settings for testing purposes
         )[design]

settings_used <- use
##' ## SET TABLES TO CREATE
scripts_used <- c("bench1a.r", "bench2a.r", "bench3a.r")[table]


library(gRbase)
library(doBy)
library(tidyverse)
library(parallel)
library(knitr)
library(DT)
library(kableExtra)
library(dplyr)
library(tidyr)
library(glue)
library(pander)
library(RhpcBLASctl)
library(gRips)

ENG = "cpp"


##' ## Dimensions etc for paper
design_paper <- list(
    NAME = "design_paper", EPS = EPS, ENG = ENG, 
    NVAR_FIX        = 100,          
    NVAR_VEC        = c(50, 100),
    PROB_VEC        = c(0.1, 0.30, 0.50, 0.70),
    NMOD            = 5,
    GRID_SCALING    = 1,
    TREE_NMOD       = 5,
    TREE_PROB       = c(0.001, 0.005, 0.010),
    TREE_USE        = 10 * c(50, 100, 200, 400, 600)
)

##' ## Dimensions etc for testing script file (small)
design_small <- list(
    NAME = "design_small", EPS = EPS, ENG = ENG, 
    NVAR_FIX         = 25,
    NVAR_VEC         = c(25, 50),
    PROB_VEC         = c(0.1, 0.30, 0.50, 0.70),
    NMOD             = 2,
    GRID_SCALING     = 4,
    TREE_NMOD        = 2,
    TREE_PROB        = c(0.001, 0.005, 0.01),    
    TREE_USE         = c(50, 100, 200, 400, 600)
)

##' ## Dimensions etc for testing script file (medium)
design_medium <- list(
    NAME = "design_medium", EPS = EPS, ENG = ENG,     
    NVAR_FIX         = 50,
    NVAR_VEC         = c(50, 100),
    PROB_VEC         = c(0.1, 0.30, 0.50, 0.70),
    NMOD             = 3,
    GRID_SCALING     = 2,
    TREE_NMOD        = 3,    
    TREE_PROB        = c(0.001, 0.005, 0.01),    
    TREE_USE         = 2 * c(50, 100, 200, 400, 600)
)

design_used <-
    list(design_small=design_small,
         design_medium=design_medium,
         design_paper=design_paper)[[use]]

NVAR_VEC     = design_used$NVAR_VEC
PROB_VEC     = design_used$PROB_VEC
NVAR_FIX     = design_used$NVAR_FIX
NMOD         = design_used$NMOD
GRID_SCALING = design_used$GRID_SCALING
TREE_NMOD    = design_used$TREE_NMOD
TREE_PROB    = design_used$TREE_PROB
TREE_UPPER   = design_used$TREE_UPPER
TREE_USE     = design_used$TREE_USE

names(NVAR_VEC)  = paste0("nvar_", NVAR_VEC)
names(PROB_VEC)  = paste0("prob_", PROB_VEC)

blas_set_num_threads(1)
omp_set_num_threads(1)

options("mc.cores"=parallel::detectCores(),
        "digits"=4, "width"=200)

##' ## output files
RDFILE1 <- "result_table1" 
RDFILE2 <- "result_table2" 
RDFILE3 <- "result_table3"

#' Create list of dataframes with data

if (!exists("data_lst")) {
    pros  <- gRips::prostate$x
    n_obs <- nrow(pros)
    set.seed(2022)
    n01   <- as.data.frame(matrix(rnorm(6033 * n_obs), nr=n_obs))
    Sn    <- cov2cor(cov(n01))
    Sp    <- cov2cor(cov(pros))
    data_lst <- list(pros=Sp, n01=Sn)
}

create_report_dir <- function(xtra="") {
    ## id <- paste0(grep("[0-9]", scripts_used),collapse="")

    bench_id <- gsub(".*([0-9]+).*$", "\\1", scripts_used)
    paste0(bench_id, collapse = "")
    st <- Sys.time() |> as.character()  |> gsub("\\..*","",x=_)
    st <- gsub(" |:","_", st)
    .REPORT_DIR <- paste0("reports/report_", Sys.info()[["sysname"]], "_", st, "_script_", bench_id)    
    if (nchar(xtra) > 0)
        .REPORT_DIR <- paste0(.REPORT_DIR, "_info_", xtra)
    ## .REPORT_DIR <- paste0("reports/report_", Sys.info()[["sysname"]], "_", gsub(" ","_", Sys.time()))
    ## .REPORT_DIR <- gsub(":","_",.REPORT_DIR)
    cat(sprintf("REPORT_DIR : %s\n", .REPORT_DIR))
    dir.create(.REPORT_DIR)
    .REPORT_DIR
}

copy_files_to_report_dir <- function(.RES_DIR, .REPORT_DIR) {
    fl <- list.files(.RES_DIR, pattern="*\\.r|*\\.md|*\\.html|*\\.RData", full.names = TRUE)
    print(fl)
    file.copy(fl, .REPORT_DIR, overwrite = TRUE)
    ## file.remove(fl)        
}

global_args <- list(eps=EPS, nobs=n_obs, maxit=50000)

options(width=140)
