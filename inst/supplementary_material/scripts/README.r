##' -----------------------------------------------------------------------
##' -----------------------------------------------------------------------
##'
##' Produce tables for the manuscript
##'
##' Højsgaard and Lauritzen (2023) On some algorithms for estimation
##' in Gaussian graphical models
##'
##' Should anything unforseen happen, please do not hesitate to contact
##' Søren Højsgaard; sorenh@math.aau.dk
##'
##' -----------------------------------------------------------------------
##' -----------------------------------------------------------------------

##' Please read the entire text here.
##'
##' To reproduce the results do the following:
##'
##' 1) Install necessary software (see below)
##'
##' 2) In the file bench_settings.r, one can specify (a) which table
##' to generate and (b) the size of the experimental design. (It takes
##' quite some time to generate the tables in the manuscript, so it
##' may be of interest to work with smaller experimental designs.)
##'
##' 3) To generate tables edit bench_settings.r as described above and
##' do source("makefile_mimic.r"). The output is a directory with the
##' results in the reports directory.
##' 


##' Install necessary software

##' a): Install R from https://www.r-project.org

##' b): For windows users: If you are on a windows platform: Install rtools
##' from https://cran.r-project.org/bin/windows/Rtools/

##' c): Prepare R for downloading packages which are hosted on
##' Bioconductor.org. To do so, run the following lines 

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("gRbase", "gRain", "gRim"), force=TRUE)

##' d): Install the necessary packages by running the following lines

pkg.list <- c(
    "gRbase",
    "doBy",
    "tidyverse",
    "knitr",
    "DT",
    "kableExtra",
    "dplyr",
    "tidyr",
    "glue",
    "pander",
    "remotes",
    "RhpcBLASctl")

sapply(pkg.list, function(pkg) {
    cat("Installing: ", pkg, "\n")
    install.packages(pkg)
})		 

##' Check that all packages are installed with:

pkg.list %in% rownames(installed.packages())

##' e): Install gRips:

remotes::install_github("hojsgaard/gRips") 


