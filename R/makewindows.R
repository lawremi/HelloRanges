### =========================================================================
### bedtools makewindows command
### -------------------------------------------------------------------------
###

bedtools_makewindows <- function(cmd = "--help") {
    do_R_call(R_bedtools_makewindows, BEDTOOLS_MAKEWINDOWS_DOC, cmd)
}

### We do not support -i and -reverse, since they seem unnecessary
### given that the window identities are implied by the list structure.

R_bedtools_makewindows <- function(b, g = NA_character_, w, s, n)
{
    stopifnot(missing(b) || isSingleString(b) || hasRanges(b),
              isSingleStringOrNA(g), !(missing(b) && missing(g)),
              missing(w) || isSingleString(w),
              missing(s) || isSingleString(s),
              missing(n) || (isSingleString(n) && missing(w) && missing(s)))

    importGenome(g)

    if (missing(b)) {
        .gr_b <- quote(as(genome, "GRanges"))
    } else {
        .gr_b <- importA(b)
    }
    
    if (!missing(s)) {
        w <- as.integer(w)
        s <- as.integer(s)
        R(ans <- slidingWindows(.gr_b, w, s))
    } else if (!missing(w)) {
        w <- as.integer(w)
        R(ans <- tile(.gr_b, width=w))
    } else {
        n <- as.integer(n)
        R(ans <- tile(.gr_b, n))
    }
    
    R(ans)
}

BEDTOOLS_MAKEWINDOWS_DOC <-
    "Usage:
       bedtools_makewindows [options]
     Options:
       -b <FILE>  BAM/BED/GFF/VCF file to tile.
       -g <path>  Genome (file or identifier) to tile.
       -w <size>  Window size.
       -s <size>  Step size for sliding windows.
       -n <count>  Number of windows, exclusive with -w.
"

do_bedtools_makewindows <- make_do(R_bedtools_makewindows)
