### =========================================================================
### utilities for defining the package
### -------------------------------------------------------------------------
###

make_do <- function(fun) {
    as.function(c(formals(fun), list(call("eval", body(fun)))))
}

