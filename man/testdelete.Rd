% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_single_edges.R
\name{testdelete}
\alias{testdelete}
\alias{testdelete.iModel}
\alias{print.testdelete}
\alias{testdelete.mModel}
\title{Test deletion of edge from an interaction model}
\usage{
testdelete(object, edge, k = 2, details = 1, ...)
}
\arguments{
\item{object}{A model; an object of class \code{iModel}.}

\item{edge}{An edge in the model; either as a right-hand sided
formula or as a vector}

\item{k}{Penalty parameter used when calculating change in AIC}

\item{details}{The amount of details to be printed; 0 surpresses
all information}

\item{\dots}{Further arguments to be passed on to the underlying
functions for testing.}
}
\value{
A list.
}
\description{
Tests if an edge can be deleted from an interaction
model.
}
\details{
If the model is decomposable and the edge is contained in
one clique only then the test is made in the marginal model
given by that clique. In that case, if the model is a
log-linear model then degrees of freedom are adjusted for
sparsity

If model is decomposable and edge is in one clique only,
then degrees of freedom are adjusted for sparsity
}
\examples{

## Discrete models
data(reinis)

## A decomposable model
mf <- ~smoke:phys:mental + smoke:systol:mental
object <- dmod(mf, data=reinis)
testdelete(object, c("phys", "mental"))
testdelete(object, c("smoke", "mental"))

## A non-decomposable model
mf <- ~smoke:phys + phys:mental + smoke:systol + systol:mental
object <- dmod(mf, data=reinis)

testdelete(object, c("phys", "mental"))

## Continuous models
data(math)

## A decomposable model
mf <- ~me:ve:al + me:al:an
object <- cmod(mf, data=math)
testdelete(object, c("ve", "al"))
testdelete(object, c("me", "al"))

## A non-decomposable model
mf <- ~me:ve + ve:al + al:an + an:me
object <- cmod(mf, data=math)
testdelete(object, c("me", "ve"))
}
\seealso{
\code{\link{testadd}}
}
\author{
Søren Højsgaard, \email{sorenh@math.aau.dk}
}
\keyword{htest}
\keyword{models}
