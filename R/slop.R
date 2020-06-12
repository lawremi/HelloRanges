### =========================================================================
### bedtools slop command
### -------------------------------------------------------------------------
###

bedtools_slop <- function(cmd = "--help") {
    do_R_call(R_bedtools_slop, BEDTOOLS_SLOP_DOC, cmd)
}

R_bedtools_slop <- function(i, b = 0, l = 0, r = 0, s = FALSE,
                            pct = FALSE, g = NULL, header = FALSE)
{
    stopifnot(isSingleString(i) || hasRanges(i),
              isSingleNumber(b), b >= 0L,
              isSingleNumber(l), l >= 0L,
              isSingleNumber(r), r >= 0L,
              xor(!(missing(l) && missing(r)), !missing(b)),
              isTRUEorFALSE(s), !(s && b),
              isTRUEorFALSE(pct),
              isGenome(g),
              isTRUEorFALSE(header))
    
    importGenome(g)
    
    i <- normA(i)
    .gr_i <- importA(i)
    .gr_i_o <- prepOverlapRanges(i, FALSE)

    if (b) {
        if (pct) {
            b <- .R(width(.gr_i_o) * b)
        }
       R(ans <- .gr_i_o + b)
    } else {
        ignore.strand <- !s
        have_l <- l != 0
        have_r <- r != 0
        if (pct) {
            if (have_l) l <- .R(width(.gr_i_o) * l)
            if (have_r) r <- .R(width(.gr_i_o) * r)
        }
        if (have_l) {
            R(ans <- resize(.gr_i_o, l + width(.gr_i_o), fix="end",
                            ignore.strand=ignore.strand))
            .gr_i_o <- quote(ans)
        }
        if (have_r) {
            R(ans <- resize(.gr_i_o, r + width(.gr_i_o),
                            ignore.strand=ignore.strand))
        }
    }

    R(ans)
}

BEDTOOLS_SLOP_DOC <-
    "Usage:
       bedtools_slop [options]
     Options:
       -i <FILE>  BAM/BED/GFF/VCF file.
       -b <size>  Increase the BED/GFF/VCF entry by the same number base pairs
                  in each direction. Integer.
       -l <size>  The number of base pairs to subtract from the start
                  coordinate. Integer.
       -r <size>  The number of base pairs to add to the end coordinate.
                  Integer.
       -s  Define -l and -r based on strand. For example. if used, -l 500
           for a negative-stranded feature, it will add 500 bp to the end
           coordinate.
    --pct  Define -l and -r as a fraction of the feature's length. E.g. if used
           on a 1000bp feature, -l 0.50, will add 500 bp \"upstream\".
           Default = false.
       -g <path>  Specify a genome file or identifier that defines the order
          and size of the sequences.
 --header  Print the header from the input file prior to results.
"

do_bedtools_slop <- make_do(R_bedtools_slop)
