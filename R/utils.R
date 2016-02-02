### =========================================================================
### internal utilities
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Language utilities
###

do_R_call <- function(R_FUN, doc, cmd) {
    opts <- docopt(doc, cmd)
    opts[grep("^-", names(opts))] <- NULL
    do.call(R_FUN, opts)
}

substituteArgs <- function(expr) {
    eval(call("substitute", expr), parent.frame(2L))
}

.R <- function(expr) {
    substituteArgs(substitute(expr))
}

R <- function(expr) {
    code <- substituteArgs(substitute(expr))
    parent.frame()$.code <- c(parent.frame()$.code, code)
    parent.frame()$.code
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Vector and Ranges utilities (push up eventually)
###

setGeneric("invert", function(x, ...) standardGeneric("invert"))

setMethod("invert", "GenomicRanges", function(x) {
              map <- c("-", "+", "*")
              strand(x) <- map[as.factor(strand(x))]
              x
          })

setMethod("merge", c("Vector", "Vector"),
          function(x, y, by = findMatches(x, y), all.x = FALSE,
                   NA.VALUE = y[NA], y.name = "y") {
              stopifnot(is(by, "Hits"),
                        isTRUEorFALSE(all.x),
                        !all.x || !missing(ifnotfound))
              ans <- x[queryHits(by)]
              mcols(ans)$y <- y[subjectHits(by)]
              if (all.x) {
                  only_x <- rep(TRUE, queryLength(by))
                  only_x[queryHits(by)] <- FALSE
                  ans_only_x <- x[only_x]
                  mcols(ans_only_x)$y <- ifnotfound
                  ans <- c(ans, ans_only_x)
              }
              ans
          })

