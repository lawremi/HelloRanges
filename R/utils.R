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

substituteArgs <- function(expr, env = parent.frame(2L)) {
    eval(call("substitute", expr), env)
}

.R <- function(expr) {
    substituteArgs(substitute(expr))
}

R <- function(expr, env = parent.frame()) {
    code <- substituteArgs(substitute(expr), env)
    env$.code <- as.call(c(quote(`{`), as.list(env$.code[-1L]), code))
    env$.code
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Genome file handling
###

setClass("GenomeFile", contains="RTLFile")

GenomeFile <- function(resource) {
    new("GenomeFile", resource=resource)
}

setMethod("import", "GenomeFile", function (con, format, text, ...) {
    if (!missing(format))
        checkArgFormat(con, format)
    df <- read.table(con, sep="\t", colClasses=c("character", "integer"),
                     col.names=c("seqnames", "seqlengths"))
    genome <- file_path_sans_ext(basename(path(con)))
    with(df, Seqinfo(seqnames, seqlengths, genome=genome))
})

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

## Once defined methods on 'merge' but generating a Pairs is not merging

setGeneric("pair", function(x, y, ...) standardGeneric("pair"))

setMethod("pair", c("Vector", "Vector"),
          function(x, y, by = findMatches(x, y), all.x = FALSE,
                   NA.VALUE = y[NA])
          {
              stopifnot(is(by, "Hits"),
                        isTRUEorFALSE(all.x))
              ans <- Pairs(x, y, hits=by)
              if (all.x) {
                  only_x <- rep(TRUE, queryLength(by))
                  only_x[queryHits(by)] <- FALSE
                  ans_only_x <- Pairs(x[only_x], rep(NA.VALUE, sum(only_x)))
                  ans <- c(ans, ans_only_x)
                  ans <- ans[order(c(queryHits(by), which(only_x)))]
              }
              ans
          })

setMethods("pair",
           list(c("GenomicRanges", "GenomicRanges"),
                c("GAlignments", "GenomicRanges"),
                c("SummarizedExperiment", "GenomicRanges")),
           function(x, y, by = findMatches(x, y), all.x = FALSE,
                    NA.VALUE = NAGRanges(y))
           {
               stopifnot(is(NA.VALUE, "GenomicRanges"))
               seqlevels(x) <- union(seqlevels(NA.VALUE), seqlevels(x))
               seqlevels(y) <- union(seqlevels(NA.VALUE), seqlevels(y))
               callNextMethod(x, y, by, all.x, NA.VALUE)
           })
