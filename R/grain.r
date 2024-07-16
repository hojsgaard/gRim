#' @title Coerce model object to Bayesian netowork
#'
#' @param x A `dmod` object
#'
#'
#' @examples
#'
#' fm <- dmod(~.^., data=lizard) |> stepwise()
#' fm |> formula()
#' bn <- as.grain(fm)
#' bn
#' 
#' @export
as.grain <- function(x){
    UseMethod("as.grain")
}

#' @export
as.grain.dModel <- function(x){
  uu <- ug(terms(x))
  out <- grain(uu, data=x$datainfo$data)  
  return(out)  
}
