## FIXME: One day also coerce to dgCMatrix; not urgent

setOldClass("cModel")
setOldClass("dModel")
setOldClass("mModel")

setAs("dModel", "matrix",      function(from) {ugList(.glist(from), result="matrix"  )})
setAs("cModel", "matrix",      function(from) {ugList(.glist(from), result="matrix"  )})
setAs("mModel", "matrix",      function(from) {ugList(.glist(from), result="matrix"  )})
setAs("dModel", "igraph",      function(from) {ugList(.glist(from), result="igraph"  )})
setAs("cModel", "igraph",      function(from) {ugList(.glist(from), result="igraph"  )})
setAs("mModel", "igraph",      function(from) {ugList(.glist(from), result="igraph"  )})

