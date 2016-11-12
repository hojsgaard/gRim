

#' Class \code{"dModel"}
#' 
#' Setting formal classes for dModel, cModel and mModel objects
#' 
#' 
#' @name dModel-class
#' @aliases dModel-class cModel-class mModel-class
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @keywords classes
#' @examples
#' 
#' showClass("dModel")
#' 
NULL





#' Internal functions in the gRim package
#' 
#' Internal functions called by other functions.
#' 
#' 
#' @aliases ghk2pmsParms_ pms2ghkParms_ normalize_ghkParms_ update_ghkParms_
#' getmi dModel cModel mModel addEdge_glist addTerm_glist dropEdge_glist
#' dropTerm_glist logLik.iModel summary.iModel print.iModelsummary
#' update.iModel extractAIC.iModel print.iModel plot.mModel plot.iModel
#' iplot.iModel iplot.mModel fit.cModel fit.dModel fit.mModel logLik.mModel
#' update.mModel formula.iModel terms.iModel CItest.data.frame
#' compareModels.iModel print.compareiModels CGstats2mmodParms moment2pmsParms
#' compareGC compareGCpairs pFormula print.MIfit print.MIparms
#' weakMarginalModel weakMarginalData
#' @keywords internal






#' Stepwise model selection in (graphical) interaction models
#' 
#' Stepwise model selection in (graphical) interaction models
#' 
#' 
#' @aliases stepwise.iModel backward forward
#' @param object An \code{iModel} model object
#' @param criterion Either \code{"aic"} or \code{"test"} (for significance
#' test)
#' @param alpha Critical value for deeming an edge to be significant/
#' insignificant. When \code{criterion="aic"}, \code{alpha} defaults to 0; when
#' \code{criterion="test"}, \code{alpha} defaults to 0.05.
#' @param type Type of models to search. Either \code{"decomposable"} or
#' \code{"unrestricted"}. If \code{type="decomposable"} and the initial model
#' is decompsable, then the search is among decomposable models only.
#' @param search Either \code{'all'} (greedy) or \code{'headlong'} (search
#' edges randomly; stop when an improvement has been found).
#' @param steps Maximum number of steps.
#' @param k Penalty term when \code{criterion="aic"}. Only k=2 gives genuine
#' AIC.
#' @param fixinMAT Matrix (p x 2) of edges. If those edges are in the model,
#' they are not considered for removal.
#' @param fixoutMAT Matrix (p x 2) of edges. If those edges are not in the
#' model, they are not considered for addition.
#' @param direction Direction for model search. Either \code{"backward"} or
#' \code{"forward"}.
#' @param details Controls the level of printing on the screen.
#' @param trace For debugging only.
#' @param \dots Further arguments to be passed on to \code{testdelete} (for
#' \code{testInEdges}) and \code{testadd} (for \code{testOutEdges}).
#' @return An \code{iModel} model object.
#' @author S<f8>ren H<f8>jsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{cmod}} \code{\link{dmod}} \code{\link{mmod}}
#' \code{\link{testInEdges}} \code{\link{testOutEdges}}
#' @keywords models
#' @examples
#' 
#' data(reinis)
#' ## The saturated model
#' m1 <- dmod(~.^., data=reinis)
#' m2 <- stepwise(m1)
#' m2
#' 
#' 






#' Return the dimension of a log-linear model
#' 
#' Return the dimension of a log-linear model given by the generating class
#' 'glist'. If the model is decomposable and adjusted dimension can be found.
#' 
#' \code{glist} can be either a list of vectors with variable names or a list
#' of vectors of variable indices.
#' 
#' \code{tableinfo} can be one of three different things.
#' 
#' 1) A contingency table (a \code{table}).
#' 
#' 2) A list with the names of the variables and their levels (such as one
#' would get if calling \code{dimnames} on a \code{table}).
#' 
#' 3) A vector with the levels. If \code{glist} is a list of vectors with
#' variable names, then the entries of the vector \code{tableinfo} must be
#' named.
#' 
#' If the model is decomposable it \code{loglinDecDim} is to be preferred over
#' \code{loglinGenDim} as the former is much faster.
#' 
#' Setting \code{adjust=TRUE} will force \code{loglinDecDim} to calculated a
#' dimension which is adjusted for sparsity of data. For this to work,
#' \code{tableinfo} *MUST* be a table.
#' 
#' @aliases loglinGenDim loglinDecDim
#' @param glist Generating class (a list) for a log-linear model. See 'details'
#' below.
#' @param tableinfo Specification of the levels of the variables. See 'details'
#' below.
#' @param adjust Should model dimension be adjusted for sparsity of data (only
#' available for decomposable models)
#' @return A numeric.
#' @author S<f8>ren H<f8>jsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{dmod}}, \code{\link{glm}}, \code{\link{loglm}}
#' @keywords models
#' @examples
#' 
#' 
#' ## glist contains variable names and tableinfo is a named vector:
#' loglinGenDim(list(c("a","b"),c("b","c")), c(a=4,b=7,c=6))
#' 
#' ## glist contains variable names and tableinfo is not named:
#' loglinGenDim(list(c(1,2),c(2,3)), c(4,7,6))
#' 
#' ## For decomposable models:
#' loglinDecDim(list(c("a","b"),c("b","c")), c(a=4,b=7,c=6),adjust=FALSE)
#' 






#' Test edges in graphical models with p-value/AIC value
#' 
#' Test edges in graphical models with p-value/AIC value. The models must
#' \code{iModel}s.
#' 
#' 
#' @aliases testInEdges testOutEdges testEdges testEdges.iModel
#' @param object An \code{iModel} model object
#' @param edgeMAT A p * 2 matrix with edges
#' @param criterion Either \code{"aic"} or \code{"test"} (for significance
#' test)
#' @param k Penalty term when \code{criterion="aic"}. Only k=2 gives genuine
#' AIC.
#' @param alpha Critical value for deeming an edge to be significant/
#' insignificant. When \code{criterion="aic"}, \code{alpha} defaults to 0; when
#' \code{criterion="test"}, \code{alpha} defaults to 0.05.
#' @param headlong If TRUE then testing will stop once a model improvement has
#' been found.
#' @param details Controls the level of printing on the screen.
#' @param \dots Further arguments to be passed on to \code{testdelete} (for
#' \code{testInEdges}) and \code{testadd} (for \code{testOutEdges}).
#' @return A matrix.
#' @author S<f8>ren H<f8>jsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{getEdges}}, \code{\link{testadd}},
#' \code{\link{testdelete}}
#' @keywords models htest
#' @examples
#' 
#' data(math)
#' cm1 <- cmod(~me:ve+ve:al+al:an, data=math)
#' testInEdges(cm1, getEdges(cm1$glist))
#' testOutEdges(cm1, getEdges(cm1$glist, ingraph=FALSE))
#' 




