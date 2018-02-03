\name{loglinDim}
\alias{loglinGenDim}
\alias{loglinDecDim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Return the dimension of a log-linear model}
\description{
  Return the dimension of a log-linear model given by the generating
  class 'glist'. If the model is decomposable and adjusted dimension can
  be found.
}
\usage{
loglinGenDim(glist, tableinfo)
loglinDecDim(glist, tableinfo, adjust=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{glist}{Generating class (a list) for a log-linear model. See
    'details' below.}
  
  \item{tableinfo}{Specification of the levels of the variables. See
    'details' below.}
  
  \item{adjust}{Should model dimension be adjusted for sparsity of data
    (only available for decomposable models)}
}
\details{
  \code{glist} can be either a list of vectors with variable names or a
  list of vectors of variable indices.
  
  \code{tableinfo} can be one of three different things.

  1) A contingency table (a \code{table}).

  2) A list with the names of the variables and their levels (such as
  one would get if calling \code{dimnames} on a \code{table}).

  3) A vector with the levels. If \code{glist} is a list of vectors with
  variable names, then the entries of the vector \code{tableinfo} must
  be named.

  If the model is decomposable it \code{loglinDecDim} is to be preferred
  over \code{loglinGenDim} as the former is much faster.

  Setting \code{adjust=TRUE} will force \code{loglinDecDim} to
  calculated a dimension which is adjusted for sparsity of data. For
  this to work, \code{tableinfo} *MUST* be a table.
}
\value{
  A numeric.
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Søren Højsgaard, \email{sorenh@math.aau.dk}
}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{
  \code{\link{dmod}},
  \code{\link{glm}},
  \code{\link{loglm}}
}
\examples{

## glist contains variable names and tableinfo is a named vector:
loglinGenDim(list(c("a","b"),c("b","c")), c(a=4,b=7,c=6))

## glist contains variable names and tableinfo is not named:
loglinGenDim(list(c(1,2),c(2,3)), c(4,7,6))

## For decomposable models:
loglinDecDim(list(c("a","b"),c("b","c")), c(a=4,b=7,c=6),adjust=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
