\name{ciTest}
\alias{ciTest}
\alias{ciTest.data.frame}
\alias{ciTest.table}
\alias{ciTest.list}

\alias{print.citest}
\alias{summary.citest}

%\alias{summary.citest}
%\alias{print.citest}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generic function for conditional independence test
}
\description{
  Generic function for conditional independence test. Specializes to
  specific types of data. 
}
\usage{
ciTest(x, set=NULL, ...)
\method{ciTest}{data.frame}(x, set=NULL, ...)
\method{ciTest}{table}(x, set=NULL, ...)
\method{ciTest}{list}(x, set=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object for which a test for conditional independence is to be
    made. See 'details' for valid types of \code{x}.
  }
  \item{set}{
    A specification of the test to be made. The tests are of the form u
    and v are independent condionally on S where u and v are variables
    and S is a set of variables. See 'details' for details about
    specification of \code{set}.
  }
  
  \item{\dots}{
    Additional arguments to be passed on to other methods.
}
}
\details{
  \code{x} can be 1) a table, 2) a dataframe whose columns are
  numerics and factors or 3) a list with components \code{cov} and
  \code{n.obs}.
  
  \code{set} can be 1) a vector or 2) a right-hand sided formula in
  which variables are separated by '+'. In either case, it is tested if
  the first two variables in the \code{set} are conditionally
  independent given the remaining variables in \code{set}. 
  (Notice an abuse of the '+' operator in the right-hand sided formula:
  The order of the variables does matter.)
  
}
\value{
  An object of class 'citest' (which is a list).
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Søren Højsgaard, \email{sorenh@math.aau.dk}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%  \code{\link{ciTest}}
  
  \code{\link{ciTest.table}}
  \code{\link{ciTest_table}}
  
  \code{\link{ciTest.data.frame}}
  \code{\link{ciTest_df}}
  
  \code{\link{ciTest.list}}
  \code{\link{ciTest_mvn}}

  \code{\link{chisq.test}}  
}
\examples{

## contingency table:
data(reinis)
## dataframe with only numeric variables:
data(carcass)
## dataframe with numeric variables and factors:
data(milkcomp1)

ciTest(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
ciTest(reinis, set=~smo+phy+sys)
ciTest(milkcomp1, set=~tre+fat+pro)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
