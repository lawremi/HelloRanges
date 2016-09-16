### =========================================================================
### bedtools nuc command
### -------------------------------------------------------------------------
###

bedtools_nuc <- function(cmd = "--help") {
    do_R_call(R_bedtools_nuc, BEDTOOLS_NUC_DOC, cmd)
}

R_bedtools_nuc <- function(fi, bed, s = FALSE, pattern = NULL,
                           fullHeader = FALSE)
{
    stopifnot(isSingleString(fi) || is(fi, "XStringSet"),
              isSingleString(bed) || hasRanges(bed),
              isTRUEorFALSE(s),
              is.null(pattern) || isSingleString(pattern),
              isTRUEorFALSE(fullHeader))
    
    if (identical(fi, "stdin"))
        stop("FASTA reading does not support 'stdin'")
    bed <- normA(bed)

### NOTE: we assume DNA here; might not be true
    R(seq <- readDNAStringSet(fi))
    if (!fullHeader) {
        R(names(seq) <- sub(" .*", "", names(seq)))
    }
    R(genome <- Seqinfo(genome=NA_character_))
    .gr_bed <- importA(bed)
    
    .gr_bed_o <- prepOverlapRanges(bed, FALSE)
    if (!s) {
        .gr_bed_o <- .R(unstrand(.gr_bed_o))
    }
    
    R(ans <- seq[.gr_bed_o])

    R(GC <- as.vector(letterFrequency(ans, "GC", as.prob=TRUE)))
    R(AT <- 1 - GC)
    R(baseCounts <- alphabetFrequency(ans, baseOnly=TRUE))

    R(mcols(ans) <- DataFrame(AT, GC, baseCounts))

    if (!is.null(pattern)) {
        R(mcols(ans)$npattern <- vcountPattern(pattern, ans))
    }
    
    R(ans)
}

BEDTOOLS_NUC_DOC <-
    "Usage:
       bedtools_nuc [options]
     Options:
       --fi <FILE>  FASTA file.
      --bed <FILE>  BED/GFF/BAM/VCF file as query.
         -s  Force strandedness. If the feature occupies the antisense strand,
             the sequence will be reverse complemented. [default: FALSE].
  --pattern <PAT>  Report the number of times a user-defined sequence
             is observed (case-sensitive).
--fullHeader Use full fasta header."

do_bedtools_nuc <- make_do(R_bedtools_nuc)
