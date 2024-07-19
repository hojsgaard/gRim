rnw2rmd <- function(infile, outfile=NULL){

    x <- readLines(infile)
    if (is.null(outfile))
        outfile <- gsub(".([r|R])md", "2.\\1md", infile)
    
    x <- gsub("\\\\verb\\|([^|]+)\\|", "`\\1`", x)
    x <- gsub("\\\\code\\{([^|]+)\\}", "`\\1`", x)
    ## x <- gsub("\\\\comi\\{([^|]+)\\}", "`\\1`", x)
    x <- gsub("\\\\comi\\{(.+?)\\}", "`\\1`", x)

    x <- gsub("%def", "", x)
    x
    writeLines(x, outfile)
}

## mdsr::Rnw2Rmd <- function (path, new_path = NULL) 
## {
##     if (is.null(new_path)) {
##         new_path <- path %>% gsub(".Rnw", ".Rmd", .) %>% gsub(".tex", 
##             ".Rmd", .) %>% fs::fs_path()
##     }
##     x <- readLines(path)
##     x <- gsub("(<<)(.*)(>>=)", "```{r \\2}", x)
##     x <- gsub("^@", "```", x)
##     x <- gsub("(\\\\Sexpr\\{)([^\\}]+)(\\})", "`r \\2`", x)
##     x <- gsub("(\\\\chapter\\{)([^\\}]+)(\\})", "# \\2", x)
##     x <- gsub("(\\\\section\\{)([^\\}]+)(\\})", "## \\2", x)
##     x <- gsub("(\\\\subsection\\{)([^\\}]+)(\\})", "### \\2", 
##         x)
##     x <- gsub("(\\\\subsubsection\\{)([^\\}]+)(\\})", "#### \\2", 
##         x)
##     x <- gsub("(\\\\citep\\{)([^\\}]+)(\\})", "[@\\2]", x)
##     x <- gsub("(\\\\cite\\{)([^\\}]+)(\\})", "@\\2", x)
##     x <- gsub("(\\\\ref\\{)([^\\}]+)(\\})", "\\\\@ref(\\2)", 
##         x)
##     x <- gsub("(\\\\label\\{)([^\\}]+)(\\})", "{#\\2}", x)
##     x <- gsub("(\\\\index\\{)([^\\}]+)(\\})(\\{)([^\\}]+)(\\})\\%", 
##         "\\\\index{\\2}{\\5}", x)
##     x <- gsub("\\\\item", "- ", x)
##     x <- gsub("(\\\\emph\\{)([^\\}]+)(\\})", "*\\2*", x)
##     x <- gsub("(\\\\textit\\{)([^\\}]+)(\\})", "*\\2*", x)
##     x <- gsub("(\\\\textbf\\{)([^\\}]+)(\\})", "**\\2**", x)
##     x <- gsub("(\\\\href\\{)([^\\}]+)(\\})(\\{)([^\\}]+)(\\})", 
##         "[\\5](\\2)", x)
##     x <- gsub("(\\\\url\\{)([^\\}]+)(\\})", "(\\2)", x)
##     x <- gsub("{\\\\tt ([a-zA-Z0-9. _()=]*)} ", "`\\1` ", x, 
##         perl = TRUE)
##     writeLines(x, new_path)
## }
