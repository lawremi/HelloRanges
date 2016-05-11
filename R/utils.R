### =========================================================================
### internal utilities
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Language utilities
###

do_R_call <- function(R_FUN, doc, cmd) {
    opts <- docopt(doc, preprocessCmd(cmd))
    opts[grep("^-", names(opts))] <- NULL
    opts <- Filter(Negate(is.null), opts)
    opts <- coerceOpts(R_FUN, opts)
    do.call(R_FUN, opts)
}

preprocessCmd <- function(cmd) {
    ## multi-character short options conflict with docopt 'stacking' feature
    gsub(" -([A-z]{2,})", " --\\1", cmd)
}

coerceOpts <- function(R_FUN, opts) {
    argclass <- vapply(formals(R_FUN)[names(opts)], class, character(1L))
    coerce <- !(argclass %in% c("NULL", "call", "name"))
    opts[coerce] <- mapply(as, opts[coerce], argclass[coerce], SIMPLIFY=FALSE)
    opts
}

substituteArgs <- function(expr) {
    eval(call("substitute", expr), parent.frame(2L))
}

.R <- function(expr) {
    substituteArgs(substitute(expr))
}

R <- function(expr) {
    code <- substituteArgs(substitute(expr))
    env <- parent.frame()
    env$.code <- as.call(c(quote(`{`), as.list(parent.frame()$.code[-1L]),
                           code))
    env$.code
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### NAGRanges: represents an NA for left outer join
###

NAGRanges <- function(x) {
    na <- GRanges(".", IRanges(0L, -1L))
    mcols(na) <- DataFrame(lapply(mcols(x), `[`, NA_integer_))
    na
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
                        !all.x || !missing(NA.VALUE))
              ans <- x[queryHits(by)]
              mcols(ans)[[y.name]] <- y[subjectHits(by)]
              if (all.x) {
                  only_x <- rep(TRUE, queryLength(by))
                  only_x[queryHits(by)] <- FALSE
                  ans_only_x <- x[only_x]
                  mcols(ans_only_x)[[y.name]] <- NA.VALUE
                  ans <- c(ans, ans_only_x)
                  ans <- ans[order(c(queryHits(by), which(only_x)))]
              }
              ans
          })

