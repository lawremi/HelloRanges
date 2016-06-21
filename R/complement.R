### =========================================================================
### bedtools complement command
### -------------------------------------------------------------------------
###

bedtools_complement <- function(cmd = "--help") {
    do_R_call(R_bedtools_complement, BEDTOOLS_COMPLEMENT_DOC, cmd)
}

R_bedtools_complement <- function(i, g) {
    stopifnot(isSingleString(i) || hasRanges(i),
              isSingleString(g))

    importGenome(g)

    i <- normA(i)
    .gr_i <- importA(i)
    .gr_i_o <- prepOverlapRanges(i, FALSE)

    R(ans <- setdiff(as(seqinfo(.gr_i), "GRanges"), unstrand(.gr_i_o)))
    R(ans)
}

BEDTOOLS_COMPLEMENT_DOC <-
    "Usage:
       bedtools_complement [options]
     Options:
       -i <FILE>  BAM/BED/GFF/VCF file.
       -g <path>  Specify a genome file or identifier that defines the order and size of the sequences."

do_bedtools_complement <- make_do(R_bedtools_complement)
