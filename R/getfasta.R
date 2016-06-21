### =========================================================================
### bedtools getfasta command
### -------------------------------------------------------------------------
###

bedtools_getfasta <- function(cmd = "--help") {
    do_R_call(R_bedtools_getfasta, BEDTOOLS_GETFASTA_DOC, cmd)
}

R_bedtools_getfasta <- function(fi, bed, s = FALSE, split = FALSE)
{
    stopifnot(isSingleString(fi) || is(fi, "XStringSet"),
              isSingleString(bed) || hasRanges(bed),
              isTRUEorFALSE(s),
              isTRUEorFALSE(split))

    if (identical(fi, "stdin"))
        stop("FASTA reading does not support 'stdin'")
    bed <- normA(bed)

### NOTE: we assume DNA here; might not be true
    R(seq <- readDNAStringSet(fi))
    R(genome <- Seqinfo(genome=NA_character_))
    .gr_bed <- importA(bed)
    
    .gr_bed_o <- prepOverlapRanges(bed, split)
    if (!s) {
        .gr_bed_o <- .R(unstrand(.gr_bed_o))
    }
    
    R(ans <- seq[.gr_bed_o])

    R(ans)
}

BEDTOOLS_GETFASTA_DOC <-
    "Usage:
       bedtools_getfasta [options]
     Options:
       --fi <FILE>  FASTA file.
      --bed <FILE>  BED/GFF/BAM/VCF file as query.
         -s  Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. [default: FALSE].
    --split  Given BED12 or BAM input, extract and concatenate the sequences from the \"blocks\" (e.g., exons)"

do_bedtools_getfasta <- make_do(R_bedtools_getfasta)
