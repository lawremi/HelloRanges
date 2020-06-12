### =========================================================================
### bedtools subtract command
### -------------------------------------------------------------------------
###

bedtools_subtract <- function(cmd = "--help") {
    do_R_call(R_bedtools_subtract, BEDTOOLS_SUBTRACT_DOC, cmd)
}

R_bedtools_subtract <- function(a, b,
                                f=1e-9, F=1e-9, r=FALSE, e=FALSE,
                                s=FALSE, S=FALSE, A=FALSE, N=FALSE,
                                g=NA_character_)
{
    stopifnot(isSingleString(a) || hasRanges(a),
              (is.character(b) && !anyNA(b) && length(b) >= 1L) || hasRanges(b),
              isSingleNumber(f), f > 0, f <= 1,
              isSingleNumber(F), F > 0, F <= 1,
              isTRUEorFALSE(r),
              isTRUEorFALSE(e),
              isTRUEorFALSE(s),
              isTRUEorFALSE(S), !(s && S),
              isTRUEorFALSE(A),
              isTRUEorFALSE(N), !(N && A),
              isGenome(g))

    importGenome(g)
    
    a <- normA(a)
    b <- normB(b)
    
    .gr_a <- importA(a)
    .gr_b <- importB(b)
    
    .gr_a_o <- prepOverlapRanges(a, FALSE)
    .gr_b_o <- prepOverlapRanges(b, FALSE)

    if (S) {
        .gr_b_o <- .R(invertStrand(.gr_b_o))
    }

    ignore.strand <- !(s || S)
    
    have_f <- !identical(f, formals(sys.function())$f)
    have_F <- !identical(F, formals(sys.function())$F)

    if (have_f || have_F || !(A || N)) {
        if (N) {
            R(hits <- findOverlaps(.gr_a_o, .gr_b_o,
                                   ignore.strand=ignore.strand))
            N_f <- f # defer the 'f' restriction
            f <- formals(sys.function())$f
            have_f <- FALSE
        } else {
            have_f <- .findOverlaps(.gr_a_o, .gr_b_o, ignore.strand, f, r, e,
                                    ret.pairs=FALSE)
        }
        
        if (have_f || have_F) {
            R(pairs <- Pairs(.gr_a_o, .gr_b_o, hits=hits))
            keep <- restrictByFraction(f, F, r, e, have_f, have_F,
                                       is_grl_a=FALSE, is_grl_b=FALSE,
                                       ignore.strand)
            R(hits <- hits[keep])
        }
        if (A) {
            R(ans <- .gr_a[countQueryHits(hits) == 0L])
        } else {
            R(toSubtract <- reduce(extractList(.gr_b_o, as(hits, "List")),
                                   ignore.strand=TRUE))
            if (N) {
                R(ans <- .gr_a[sum(width(toSubtract)) / width(.gr_a_o) <= N_f])
            } else {
                R(ans <- psetdiff(.gr_a_o, toSubtract, ignore.strand=TRUE))
                R(ans <- subset(ans, width > 0L))
            }
        }
    } else {
        if (identical(.gr_a_o, .gr_a)) {
            R(ans <- subsetByOverlaps(.gr_a_o, .gr_b_o, invert=TRUE,
                                      ignore.strand=ignore.strand))
        } else {
            if (ignore.strand) {
                .gr_b_o <- .R(unstrand(.gr_b_o))
            }
            R(ans <- .gr_a[.gr_a_o %outside% .gr_b_o])
        }
    }

    R(ans)
}

BEDTOOLS_SUBTRACT_DOC <-
    "Usage:
       bedtools_subtract [options]
     Options:
       -a <FILE>  BAM/BED/GFF/VCF file A. Each feature in A is compared to B in
          search of overlaps. Use 'stdin' if passing A with a UNIX pipe.
       -b <FILE1,...>  One or more BAM/BED/GFF/VCF file(s) B. Use 'stdin' if
          passing B with a UNIX pipe. -b may be followed with multiple
          databases and/or wildcard (*) character(s).
       -f <frac>  Minimum overlap required as a fraction of A [default: 1e-9].
       -F <frac>  Minimum overlap required as a fraction of B [default: 1e-9].
       -r  Require that the fraction of overlap be reciprocal for A and B.
           In other words, if -f is 0.90 and -r is used, this requires that B
           overlap at least 90% of A and that A also overlaps at least 90% of B.
       -e  Require that the minimum fraction be satisfied for A _OR_ B. In
           other words, if -e is used with -f 0.90 and -F 0.10 this requires
           that either 90% of A is covered OR 10% of B is covered. Without -e,
           both fractions would have to be satisfied.
       -s  Force strandedness. That is, only report hits in B that overlap A on
           the same strand. By default, overlaps are reported without respect
           to strand.
       -S  Require different strandedness. That is, only report hits in B that
           overlap A on the _opposite_ strand. By default, overlaps are
           reported without respect to strand.
       -A  Remove entire feature if any overlap. That is, by default, only
           subtract the portion of A that overlaps B. Here, if any overlap is
           found (or -f amount), the entire feature is removed.
       -N  Same as -A except when used with -f, the amount is the sum of all
           features (not any single feature).
       -g <path>  Specify a genome file or identifier that defines the order
          and size of the sequences."

do_bedtools_subtract <- make_do(R_bedtools_subtract)
