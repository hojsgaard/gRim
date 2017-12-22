\name{stepwise.iModel; backward; forward}
\alias{stepwise.iModel}
\alias{backward}
\alias{forward}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stepwise model selection in (graphical) interaction models}
\description{
Stepwise model selection in (graphical) interaction models
}
\usage{
\method{stepwise}{iModel}(object, 
 criterion = "aic", alpha = NULL, type ="decomposable",
 search="all", steps = 1000, k = 2,
 direction = "backward", fixinMAT=NULL, fixoutMAT=NULL,
 details = 0, trace = 2, ...)

backward(object,
 criterion = "aic", alpha = NULL, type = "decomposable",
 search="all", steps = 1000, k = 2,
 fixinMAT=NULL, details = 1, trace = 2,...)

forward(object,
 criterion = "aic", alpha = NULL, type = "decomposable",
 search="all", steps = 1000, k = 2,
 fixoutMAT=NULL,  details = 1, trace = 2,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{iModel} model object}
  
  \item{criterion}{Either \code{"aic"} or \code{"test"} (for
    significance test)}
  
  \item{alpha}{Critical value for deeming an edge to be significant/
    insignificant. When \code{criterion="aic"}, \code{alpha} defaults to
    0; when \code{criterion="test"}, \code{alpha} defaults to 0.05.}
  
  \item{type}{Type of models to search. Either \code{"decomposable"} or
    \code{"unrestricted"}. If \code{type="decomposable"} and the initial
    model is decompsable, then the search is among decomposable models
    only. }

  \item{search}{Either \code{'all'} (greedy) or \code{'headlong'}
    (search edges randomly; stop when an improvement has been found).}

  \item{steps}{Maximum number of steps. }

  \item{k}{Penalty term when \code{criterion="aic"}. Only k=2 gives
    genuine AIC.}  

  \item{fixinMAT}{Matrix (p x 2) of edges. If those edges are in
    the model, they are not considered for removal.}
  
  \item{fixoutMAT}{Matrix (p x 2) of edges. If those edges are not
    in the model, they are not considered for addition.} 

  \item{direction}{Direction for model search. Either \code{"backward"}
    or \code{"forward"}.} 
  
  \item{details}{Controls the level of printing on the screen. }
  \item{trace}{For debugging only.}
  \item{\dots}{Further arguments to be passed on to \code{testdelete}
    (for \code{testInEdges}) and \code{testadd} (for \code{testOutEdges}). }
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  An \code{iModel} model object.
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Søren Højsgaard, \email{sorenh@math.aau.dk}
}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{
  \code{\link{cmod}}
  \code{\link{dmod}}
  \code{\link{mmod}}
  \code{\link{testInEdges}}
  \code{\link{testOutEdges}}
}
\examples{
data(reinis)
## The saturated model
m1 <- dmod(~.^., data=reinis)
m2 <- stepwise(m1)
m2

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
