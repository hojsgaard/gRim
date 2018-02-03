\name{testInEdges; testOutEdges}
\alias{testInEdges}
\alias{testOutEdges}

\alias{testEdges}
\alias{testEdges.iModel}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test edges in graphical models with p-value/AIC value}
\description{
  Test edges in graphical models with p-value/AIC value. The models
  must \code{iModel}s. 
}
\usage{
testInEdges (object, edgeMAT=NULL, criterion = "aic", k = 2,
 alpha = NULL, headlong = FALSE, details = 1, ...)
testOutEdges(object, edgeMAT=NULL, criterion = "aic", k = 2,
 alpha = NULL, headlong = FALSE, details = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{iModel} model object}
  \item{edgeMAT}{A p * 2 matrix with edges}
  \item{criterion}{Either \code{"aic"} or \code{"test"} (for significance test)}
  \item{k}{Penalty term when \code{criterion="aic"}. Only k=2 gives
    genuine AIC.}  
  \item{alpha}{Critical value for deeming an edge to be significant/
    insignificant. When \code{criterion="aic"}, \code{alpha} defaults to
    0; when \code{criterion="test"}, \code{alpha} defaults to 0.05.}
  \item{headlong}{If TRUE then testing will stop once a model
    improvement has been found. }
  \item{details}{Controls the level of printing on the screen. }
  \item{\dots}{Further arguments to be passed on to \code{testdelete}
    (for \code{testInEdges}) and \code{testadd} (for \code{testOutEdges}). }
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  A matrix.
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Søren Højsgaard, \email{sorenh@math.aau.dk}
}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{
  \code{\link{getEdges}},
  \code{\link{testadd}},
  \code{\link{testdelete}}
}
\examples{
data(math)
cm1 <- cmod(~me:ve + ve:al + al:an, data=math)
testInEdges(cm1, getEdges(cm1$glist))
testOutEdges(cm1, getEdges(cm1$glist, ingraph=FALSE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{htest}% __ONLY ONE__ keyword per line
