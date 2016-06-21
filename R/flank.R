### =========================================================================
### bedtools flank command
### -------------------------------------------------------------------------
###

bedtools_flank <- function(cmd = "--help") {
    do_R_call(R_bedtools_flank, BEDTOOLS_FLANK_DOC, cmd)
}

R_bedtools_flank <- function(i, b = 0, l = 0, r = 0, s = FALSE,
                             pct = FALSE, g = NULL, header = FALSE)
{
    stopifnot(isSingleString(i) || hasRanges(i),
              isSingleNumber(b), b >= 0L,
              isSingleNumber(l), l >= 0L,
              isSingleNumber(r), r >= 0L,
              xor(!(missing(l) && missing(r)), !missing(b)),
              isTRUEorFALSE(s), !(s && b),
              isTRUEorFALSE(pct))
    
    importGenome(g)
    
    i <- normA(i)
    .gr_i <- importA(i)
    .gr_i_o <- prepOverlapRanges(i, FALSE)

    if (b) {
        l <- b
        r <- b
    }

    ignore.strand <- !s

    if (l > 0L) {
        if (pct) {
            l <- .R(width(.gr_i_o) * l)
        }
        R(left <- flank(.gr_i_o, l, ignore.strand=ignore.strand))
    }
    if (r > 0L) {
        if (pct) {
            r <- .R(width(.gr_i_o) * r)
        }
        R(right <- flank(.gr_i_o, r, start=FALSE, ignore.strand=ignore.strand))
        if (l > 0L) {
            R(ans <- zipup(Pairs(left, right)))
        } else {
            R(ans <- right)
        }
    } else {
        R(ans <- left)
    }

    R(ans)
}

BEDTOOLS_FLANK_DOC <-
    "Usage:
       bedtools_flank [options]
     Options:
       -i <FILE,...>  BAM/BED/GFF/VCF files.
       -b <size>  Increase the BED/GFF/VCF entry by the same number base pairs in each direction. Integer.
       -l <size> The number of base pairs to subtract from the start coordinate. Integer.
       -r <size> The number of base pairs to add to the end coordinate. Integer.
       -s  Define -l and -r based on strand. For example. if used, -l 500 for a negative-stranded feature, it will add 500 bp to the end coordinate.
    --pct  Define -l and -r as a fraction of the feature's length. E.g. if used on a 1000bp feature, -l 0.50, will add 500 bp \"upstream\". Default = false.
       -g <path>  Specify a genome file or identifier that defines the order and size of the sequences.
 --header  Print the header from the input file prior to results."

do_bedtools_flank <- make_do(R_bedtools_flank)
