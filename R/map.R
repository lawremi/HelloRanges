### =========================================================================
### bedtools map command
### -------------------------------------------------------------------------
###

bedtools_map <- function(cmd = "--help") {
    do_R_call(R_bedtools_map, BEDTOOLS_MAP_DOC, cmd)
}

## NOTE: we do not support the 'null' argument; rather, we just let R
## do the sensible thing given the statistic.

R_bedtools_map <- function(a, b, c="5", o="sum",
                           f=1e-9, F=1e-9, r=FALSE, e=FALSE, s=FALSE, S=FALSE,
                           header=FALSE, # ignored
                           split=FALSE, g=NA_character_)
{
    stopifnot(isSingleString(a),
              is.character(b), !anyNA(b), length(b) >= 1L,
              is.numeric(c) || isSingleString(c),
              isSingleString(o),
              isSingleNumber(f), f > 0, f <= 1,
              isSingleNumber(F), F > 0, F <= 1,
              isTRUEorFALSE(r),
              isTRUEorFALSE(e),
              isTRUEorFALSE(s),
              isTRUEorFALSE(S), !(s && S),
              isTRUEorFALSE(split),
              isSingleStringOrNA(g),
              isTRUEorFALSE(header))

    importGenome(g)
    
    a <- normA(a)
    b <- normB(b)

    co <- normCandO(b, c, o)

    .gr_a <- importA(a)
    .gr_b <- importB(b, extraCols=co$cn)
    if (isVcf(b)) {
        .gr_b <- .R(granges(expand(.gr_b)))
    }

    .gr_a_o <- prepOverlapRanges(a, split)
    .gr_b_o <- prepOverlapRanges(b, split)

    is_grl_a <- split && (isBed(a) || isBam(a))
    is_grl_b <- split && (isBed(b) || isBam(b))

    if (S) {
        .gr_b_o <- .R(invertStrand(.gr_b_o))
    }

    ignore.strand <- !(s || S)

    have_f <- .findOverlaps(pairs=FALSE, f, r, e)
    have_F <- !identical(F, formals(sys.function())$F)

    if (have_f || have_F) {
        R(pairs <- Pairs(.gr_a_o, .gr_b_o, hits=hits))
        restrictByFraction(f, F, r, e, have_f, have_F, is_grl_a, is_grl_b,
                           ignore.strand)
        R(hits <- hits[keep])
    }

    R(ans <- .gr_a)
    
    .agg <- .aggregateCall(.gr_b, quote(hits), co$exprs, drop=FALSE)
    R(mcols(ans) <- .agg)

    R(ans)
}

BEDTOOLS_MAP_DOC <-
    "Usage:
       bedtools_map [options]
     Options:
     -a <FILE>  BAM/BED/GFF/VCF file A. Each feature in A is compared to B in search of overlaps. Use 'stdin' if passing A with a UNIX pipe.
     -b <FILE1,...> One or more BAM/BED/GFF/VCF file(s) B. Use 'stdin' if passing B with a UNIX pipe. -b may be followed with multiple databases and/or wildcard (*) character(s).
     -c <col>  Specify the column from the B file to map onto intervals in A.
         [default: 5]
     -o <op>  Specify the operation that should be applied to -c.
         Valid operations:
         sum - numeric only
         count - numeric or text
         count_distinct - numeric or text
         min - numeric only
         max - numeric only
         absmin - numeric only
         absmax - numeric only
         mean - numeric only
         median - numeric only
         collapse (i.e., print a comma separated list) - numeric or text
         distinct (i.e., print a comma separated list) - numeric or text
         concat (i.e., print a comma separated list) - numeric or text
         [default: sum]
     -f <frac>  Minimum overlap required as a fraction of A. [default: 1e-9].
     -F <frac>  Minimum overlap required as a fraction of B. [default: 1e-9].
     -r  Require that the fraction of overlap be reciprocal for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
     -e  Require that the minimum fraction be satisfied for A _OR_ B. In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of B is covered. Without -e, both fractions would have to be satisfied.
     -s  Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
     -S  Require different strandedness. That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.
-header	 Print the header from the A file prior to results.
 -split	 Treat “split” BAM (i.e., having an “N” CIGAR operation) or BED12 entries as distinct BED intervals. When using -sorted, memory usage remains low even for very large files.
     -g	 Specify a genome file or identifier the defines the expected chromosome order in the input files."

do_bedtools_map <- make_do(R_bedtools_map)
