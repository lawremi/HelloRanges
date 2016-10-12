### =========================================================================
### bedtools coverage command
### -------------------------------------------------------------------------
###

bedtools_coverage <- function(cmd = "--help") {
    do_R_call(R_bedtools_coverage, BEDTOOLS_COVERAGE_DOC, cmd)
}

R_bedtools_coverage <- function(a, b, hist=FALSE, d=FALSE, counts=FALSE,
                                f=1e-9, F=1e-9, r=FALSE, e=FALSE,
                                s=FALSE, S=FALSE,
                                split=FALSE, g=NA_character_,
                                header=FALSE, #ignored
                                sortout=FALSE)
{
    stopifnot(isSingleString(a) || hasRanges(a),
              (is.character(b) && !anyNA(b) && length(b) >= 1L) ||
                  hasRanges(b),
              isTRUEorFALSE(hist),
              isTRUEorFALSE(d),
              isTRUEorFALSE(counts), !(d && counts), !(hist && (d || counts)),
              isSingleNumber(f), f > 0, f <= 1,
              isSingleNumber(F), F > 0, F <= 1,
              isTRUEorFALSE(r),
              isTRUEorFALSE(e),
              isTRUEorFALSE(s),
              isTRUEorFALSE(S), !(s && S),
              isTRUEorFALSE(split),
              isSingleStringOrNA(g),
              isTRUEorFALSE(header),
              isTRUEorFALSE(sortout))

    importGenome(g)
    
    a <- normA(a)
    b <- normB(b)
    
    .gr_a <- importA(a)
    .gr_b <- importB(b)

    .gr_a_o <- prepOverlapRanges(a)
    .gr_b_o <- prepOverlapRanges(b, split)
    if (split && isBam(b)) {
        .gr_b_o$drop.D.ranges <- TRUE
    }

    is_grl_b <- split && (isBed(b) || isBam(b))

    if (S) {
        .gr_b_o <- .R(invertStrand(.gr_b_o))
    }

    ignore.strand <- !(s || S)

    have_f <- !identical(f, formals(sys.function())$f)
    have_F <- !identical(F, formals(sys.function())$F)
    
    if (!(have_f || have_F) && (hist || d || counts)) {
        if (counts) {
            R(ans <- .gr_a)
            R(mcols(ans)$count <- countOverlaps(.gr_a_o, .gr_b_o,
                                                ignore.strand=ignore.strand))
        } else if (ignore.strand) {
### FIXME: drop unname() once [,List,GRanges is fixed
            R(cov <- unname(coverage(.gr_b_o)[.gr_a_o]))
        } else {
            R(cov_plus <- coverage(subset(.gr_b_o, strand=="+"))[
                  subset(.gr_a_o, strand=="+")])
            R(cov_minus <- coverage(subset(.gr_b_o, strand=="-"))[
                  subset(.gr_a_o, strand=="-")])
            R(cov_star <- coverage(.gr_b_o)[subset(exp, strand=="*")])
            R(cov <- unname(unsplit(List(cov_plus, cov_minus, cov_star),
                                    strand(.gr_a_o))))
        }
    } else {
        have_f <- .findOverlaps(.gr_a_o, .gr_b_o, ignore.strand, f, r, e,
                                ret.pairs=FALSE)
        R(pairs <- Pairs(.gr_a_o, .gr_b_o, hits=hits))
        if (have_f || have_F) {
            keep <- restrictByFraction(f, F, r, e, have_f, have_F,
                                       FALSE, is_grl_b, ignore.strand)
            R(hits <- hits[keep])
            if (counts) {
                R(ans <- .gr_a)
                mcols(ans)$count <- countQueryHits(hits)
            } else {
                R(pairs <- pairs[keep])
            }
        }
        if (!counts) {
            R(int <- pintersect(pairs, ignore.strand=ignore.strand))
            R(intList <- relist(int, as(hits, "List")))
            if (hist || d) {
                R(cov <- coverage(int))
            } else {
                R(ans <- .gr_a)
                R(mcols(ans) <- within(mcols(ans), {
                    count <- countQueryHits(hits)
                    covered <- sum(width(reduce(intList)))
                    fraction <- covered / width(ans)
                }))
            }
        }
    }
    
    if (hist) {
        R(tab <- t(table(cov)))
        ## addmargins() inefficient, but rowSums() returns double...
        R(tab <- cbind(tab, all=rowSums(tab)))
        R(covhist <- DataFrame(as.table(tab)))
        R(colnames(covhist) <- c("coverage", "a", "count"))
        R(len <- c(lengths(cov, use.names=FALSE), sum(lengths(cov))))
        R(covhist$len <- rep(len, each=nrow(tab)))
        R(covhist <- subset(covhist, count > 0L))
        R(covhist$fraction <- with(covhist, count / len))
        R(ans <- .gr_a)
        rm(split, a)
        R(covhistList <- split(covhist, ~a)[,-2L])
        R(mcols(ans)$coverage <- head(covhistList, -1L))
        R(metadata(ans)$coverage <- covhistList$all)
    } else if (d) {
        R(ans <- .gr_a)
        R(mcols(ans)$coverage <- cov)
    }
   
    R(ans)
}

BEDTOOLS_COVERAGE_DOC <-
    "Usage:
       bedtools_coverage [options]
     Options:
       -a <FILE>  BAM/BED/GFF/VCF file A. Each feature in A is compared to B
          in search of overlaps. Use 'stdin' if passing A with a UNIX pipe.
       -b <FILE1,...>  One or more BAM/BED/GFF/VCF file(s) B. Use 'stdin' if
          passing B with a UNIX pipe. -b may be followed with multiple
          databases and/or wildcard (*) character(s).
   --hist  Report a histogram of coverage for each feature in A as well as a
           summary histogram for _all_ features in A.
           Output (tab delimited) after each feature in A:
              1) depth
              2) # bases at depth
              3) size of A
              4) % of A at depth
       -d  Report the depth at each position in each A feature.
           Positions reported are one based. Each position and depth
           follow the complete A feature.
 --counts  Only report the count of overlaps, don't compute fraction, etc.
           Restricted by -f and -r.
       -f <frac>  Minimum overlap required as a fraction of A [default: 1e-9].
       -F <frac>  Minimum overlap required as a fraction of B [default: 1e-9].
       -r  Require that the fraction of overlap be reciprocal for A and B.
           In other words, if -f is 0.90 and -r is used, this requires that B
           overlap at least 90% of A and A also overlaps at least 90% of B.
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
       --split  Treat split BAM (i.e., having an 'N' CIGAR operation) or BED12
                entries as distinct BED intervals.
       -g <path>  Specify a genome file or identifier that defines the order
                  and size of the sequences.
       --sortout  When using multiple databases (-b), sort the output DB hits
                  for each record.
"

do_bedtools_coverage <- make_do(R_bedtools_coverage)
