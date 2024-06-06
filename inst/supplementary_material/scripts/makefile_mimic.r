args = commandArgs(trailingOnly=TRUE)
xtra <- paste0(args, collapse = "_")

if (identical(xtra, "")) xtra <- "default"

.RES_DIR     <- tempdir()
file.remove(list.files(.RES_DIR, full=T))
cat(".RES_DIR: ", .RES_DIR, "\n")

## WHICH SETTINGS TO USE

file.copy("bench_settings.r", .RES_DIR)
list.files(.RES_DIR)

local_bench_settings <- file.path(.RES_DIR, "bench_settings.r")
source(local_bench_settings)

.REPORT_DIR <- create_report_dir(xtra)
list.files(.REPORT_DIR)

file.copy(scripts_used, .RES_DIR)
local_scripts_used   <- file.path(.RES_DIR, scripts_used) 

.t.start <- Sys.time()
cat(sprintf(".t.start: %s\n",.t.start))

save(.t.start,      file=file.path(.RES_DIR, "t.start.RData"))
save(settings_used, scripts_used,
     file=file.path(.RES_DIR, "settings_used.RData"))

parallel::mclapply(local_scripts_used,
                   function(f) {
                       source(local_bench_settings)    
                       source(f)
                   })

list.files(.RES_DIR)
copy_files_to_report_dir(.RES_DIR, .REPORT_DIR)
list.files(.RES_DIR)  
list.files(.REPORT_DIR)

if (file.exists(file.path(.REPORT_DIR, "t.start.RData"))) {
    load(file.path(.REPORT_DIR, "t.start.RData"))
    .t.end <- Sys.time()
    dtime <- difftime(.t.end, .t.start, units="mins")
    cat(paste0("minutes: ",  round(dtime,2)), file=file.path(.REPORT_DIR, "computing-time.txt"))
    ## file.remove(file.path(.REPORT_DIR, "t.start.RData"))
}


rmarkdown::render('bench_report.rmd', output_dir = .REPORT_DIR)

## move:
load(file=file.path(.REPORT_DIR, "t.start.RData"))

## create html-file with overview
rmarkdown::render("overview_bench.rmd")

.t.end <- Sys.time()
.dtime <- as.numeric(difftime(.t.end, .t.start, units="secs"))

cat(sprintf(".t.start : %s\n", .t.start))
cat(sprintf(".t.end   : %s\n", .t.end))
cat(sprintf(".dtime   : %s\n", .dtime))

cat(sprintf("TIME FOR SCRIPT %f secs %f mins, %f hrs\n",
            .dtime, .dtime/60, .dtime/3600))


## Otherwise github might choke
if (file.exists(".RData"))
    file.remove(".RData") 
