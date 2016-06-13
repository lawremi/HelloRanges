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

.R <- function(expr, env = parent.frame()) {
    substituteArgs(substitute(expr), env)
}

R <- function(expr, env = parent.frame()) {
    pushR(substituteArgs(substitute(expr), env), env)
}

pushR <- function(code = getR(parent.frame()), env = parent.frame()) {
    env$.code <- as.call(c(quote(`{`), getR(env), code))
    env$.code
}

getR <- function(env) {
    as.list(env$.code[-1L])
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
    df <- read.table(path(con), sep="\t", colClasses=c("character", "integer"),
                     col.names=c("seqnames", "seqlengths"))
    genome <- file_path_sans_ext(basename(path(con)))
    with(df, Seqinfo(seqnames, seqlengths, genome=genome))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### NAGRanges: represents an NA for left outer join
###

makeNAColumn <- function(x) {
### FIXME: extractROWS() should support NA subscripts, but oh well
    if (is(x, "Rle"))
        x <- decode(x)
    ans <- x[NA_integer_]
    if (is(x, "Rle"))
        ans <- Rle(ans)
    ans
}

NAGRanges <- function(x) {
    na <- GRanges(".", IRanges(0L, -1L))
    mcols(na) <- if (length(mcols(x)) > 0L)
                     DataFrame(lapply(mcols(x), makeNAColumn))
                 else S4Vectors:::make_zero_col_DataFrame(1L)
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mode/antimode
###

## Child and parent partition the same space, but the child
## breakpoints are nested within the parent breakpoints, i.e., for
## every parent breakpoint, there is a matching child breakpoint.
toparent <- function(child, parent) {
    PartitioningByEnd(findInterval(end(PartitioningByEnd(parent)),
                                   end(PartitioningByEnd(child))))
}

setGeneric("distmode", function(x, anti=FALSE) standardGeneric("distmode"),
           signature = "x")

setMethod("distmode", "CompressedAtomicList", function(x, anti=FALSE) {
    g <- IRanges:::subgrouping(x)
    pg <- PartitioningByEnd(g)
    counts <- relist(lengths(g), toparent(g, x))
    m <- if (anti) which.min(counts, global=TRUE)
         else which.max(counts, global=TRUE)
    unlist(x, use.names=FALSE)[unlist(g)[end(pg)[m]]]
})

setMethod("distmode", "SimpleList", function(x, anti=FALSE) {
    unlist(endoapply(x, distmode, anti=anti))
})

setMethod("distmode", "vector", function(x, anti=FALSE) {
    sm <- selfmatch(x)
    counts <- tabulate(sm)
    pos <- if (anti) which.min(counts) else which.max(counts)
    x[pos]
})

setMethod("distmode", "factor", function(x, anti=FALSE) {
    counts <- table(x)
    names(if (anti) which.min(counts) else which.max(counts))
})

