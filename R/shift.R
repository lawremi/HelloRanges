### =========================================================================
### bedtools shift command
### -------------------------------------------------------------------------
###

bedtools_shift <- function(cmd = "--help") {
    do_R_call(R_bedtools_shift, BEDTOOLS_SHIFT_DOC, cmd)
}

R_bedtools_shift <- function(i, s = 0, m = 0, p = 0,
                             pct = FALSE, g = NULL, header = FALSE)
{
    stopifnot(isSingleString(i) || hasRanges(i),
              isSingleNumber(s),
              isSingleNumber(m),
              isSingleNumber(p),
              xor(!(missing(m) && missing(p)), !missing(s)),
              isTRUEorFALSE(pct))
    
    importGenome(g)
    
    i <- normA(i)
    .gr_i <- importA(i)
    .gr_i_o <- prepOverlapRanges(i, FALSE)

    if (s != 0) {
        if (pct) {
            s <- .R(width(.gr_i_o) * s)
        }
        R(ans <- shift(.gr_i_o, s))
    } else {
        R(ans <- .gr_i_o)
    }
    if (p != 0) {
        .plus <- .R(strand(ans) == "+")
        if (pct) {
            p <- .R(width(ans) * p)
        }
        R(ans[.plus] <- shift(ans[.plus], p))
    }
    if (m != 0) {
        .minus <- .R(strand(ans) == "-")
        if (pct) {
            m <- .R(width(ans) * m)
        }
        R(ans[.minus] <- shift(ans[.minus], m))
    }

    R(ans)
}

BEDTOOLS_SHIFT_DOC <-
    "Usage:
       bedtools_shift [options]
     Options:
       -i <FILE,...>  BAM/BED/GFF/VCF files.
       -s <bp>  Shift the BED/GFF/VCF entry -s base pairs. Integer or Float (e.g. 0.1) if used with -pct.
       -m <bp>  Shift entries on the - strand -m base pairs. Integer or Float (e.g. 0.1) if used with -pct.
       -p <bp>  Shift entries on the + strand -p base pairs. Integer or Float (e.g. 0.1) if used with -pct.
    --pct  Define -l and -r as a fraction of the feature's length. E.g. if used on a 1000bp feature, -l 0.50, will add 500 bp \"upstream\". Default = false.
       -g <path>  Specify a genome file or identifier that defines the order and size of the sequences.
 --header  Print the header from the input file prior to results."

do_bedtools_shift <- make_do(R_bedtools_shift)
