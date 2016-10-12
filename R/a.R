### =========================================================================
### utilities for defining the package
### -------------------------------------------------------------------------
###

make_do <- function(fun) {
    args <- formals(fun)
    args[] <- lapply(names(args), as.name)
    as.function(c(formals(fun),
                  call("eval", as.call(c(substitute(fun), args)))))
}

globalVariables(c("count", # output subset() expression
                  "len" # output within() expression
                  ))

