### =========================================================================
### bedtools jaccard command
### -------------------------------------------------------------------------
###
### Sort of a special case tool but nicely demonstrates set operations.
###

bedtools_jaccard <- function(cmd = "--help") {
    do_R_call(R_bedtools_jaccard, BEDTOOLS_JACCARD_DOC, cmd)
}

R_bedtools_jaccard <- function(a, b,
                               f=1e-9, F=1e-9, r=FALSE, e=FALSE,
                               s=FALSE, S=FALSE,
                               split=FALSE)
{
    stopifnot(isSingleString(a),
              is.character(b), !anyNA(b), length(b) >= 1L,
              isSingleNumber(f), f > 0, f <= 1,
              isSingleNumber(F), F > 0, F <= 1,
              isTRUEorFALSE(r),
              isTRUEorFALSE(e),
              isTRUEorFALSE(s),
              isTRUEorFALSE(S), !(s && S),
              isTRUEorFALSE(split))

    R(genome <- NA_character_)
    
    a <- normA(a)
    b <- normB(b)
    
    .gr_a <- importA(a)
    .gr_b <- importB(b)
    
    .gr_a_o <- prepOverlapRanges(a, split)
    .gr_b_o <- prepOverlapRanges(b, split)

    is_grl_a <- split && (isBed(a) || isBam(a))
    is_grl_b <- split && (isBed(b) || isBam(b))

    if (S) {
        .gr_b_o <- .R(invertStrand(.gr_b_o))
    }

    ignore.strand <- !(s || S)

    have_f <- !identical(f, formals(sys.function())$f)
    have_F <- !identical(F, formals(sys.function())$F)
    
    if (have_f || have_F) {
        have_f <- .findOverlaps(pairs=TRUE, f, r, e)
        if (have_f || have_F) {
            restrictByFraction(f, F, r, e, have_f, have_F, is_grl_a, is_grl_b,
                               ignore.strand, strict.strand=TRUE)
        }
        .olap <- .R(olap[keep])
        if (is_grl_a || is_grl_b)
            .olap <- .R(unlist(.olap))
        R(intersects <- reduce(.olap, ignore.strand=ignore.strand))
    }

    if (is_grl_a) {
        .gr_a_o <- .R(unlist(.gr_a_o))
    }
    
    if (is_grl_b) {
        .gr_b_o <- .R(unlist(.gr_b_o))
    }

    if (!fracRestriction) {
        R(intersects <- intersect(.gr_a_o, .gr_b_o,
                                  ignore.strand=ignore.strand))
    }
    
    R(intersection <- sum(width(intersects)))    
    R(union <- sum(width(union(.gr_a_o, .gr_b_o, ignore.strand=ignore.strand))))
    
    R(ans <- DataFrame(intersection, union, jaccard=intersection/union,
                       n_intersections=length(intersects)))
    R(ans)
}

BEDTOOLS_JACCARD_DOC <-
    "Usage:
       bedtools_jaccard [options]
     Options:
       -a <FILE>  BAM/BED/GFF/VCF file A. Each feature in A is compared to B in search of overlaps. Use 'stdin' if passing A with a UNIX pipe.
       -b <FILE1,...> One or more BAM/BED/GFF/VCF file(s) B. Use 'stdin' if passing B with a UNIX pipe. -b may be followed with multiple databases and/or wildcard (*) character(s).
       -f <frac>  Minimum overlap required as a fraction of A [default: 1e-9].
       -F <frac>  Minimum overlap required as a fraction of B [default: 1e-9].
       -r  Require that the fraction of overlap be reciprocal for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
       -e  Require that the minimum fraction be satisfied for A _OR_ B. In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of B is covered. Without -e, both fractions would have to be satisfied.
       -s  Force strandedness. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
       -S  Require different strandedness. That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.
       --split  Treat split BAM (i.e., having an 'N' CIGAR operation) or BED12 entries as distinct BED intervals.
"

do_bedtools_jaccard <- make_do(R_bedtools_jaccard)
