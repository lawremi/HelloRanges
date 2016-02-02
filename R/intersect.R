### =========================================================================
### bedtools intersect command
### -------------------------------------------------------------------------
###

bedtools_intersect <- function(cmd = "--help") {
    do_R_call(R_bedtools_intersect, BEDTOOLS_INTERSECT_DOC, cmd)
}

stdinFile <- function() {
    message("Assuming BED format for 'stdin'; modify if otherwise")
    .R(BEDFile("stdin"))
}

normA <- function(a) {
    if (a == "stdin")
        stdinFile()
    else a
}

normB <- function(b) {
    normBToken <- function(bi) {
        if (bi == "stdin")
            stdinFile()
        else if (grepl("*", bi, fixed=TRUE))
            .R(Sys.glob(bi))
        else bi
    }
    if (length(b) == 1L)
        normToken(b)
    else as.call(c(quote(c), lapply(p, normBToken)))
}

R_bedtools_intersect <- function(a, b, ubam=FALSE, bed=FALSE,
                                 wa=FALSE, wb=FALSE, loj=FALSE, wo=FALSE,
                                 wao=FALSE, u=FALSE, c=FALSE, v=FALSE,
                                 f=1e-9, F=FALSE, r=FALSE, e=FALSE, s=FALSE,
                                 S=FALSE, split=FALSE, g=FALSE, header=FALSE,
                                 names=as.character(seq_along(b)),
                                 filenames=FALSE, sortout=FALSE)
{
    stopifnot(isSingleString(a),
              is.character(b), !anyNA(b), length(b) >= 1L,
              isTRUEorFALSE(ubam),
              isTRUEorFALSE(bed),
              isTRUEorFALSE(wa),
              isTRUEorFALSE(wb),
              isTRUEorFALSE(loj),
              isTRUEorFALSE(wo),
              isTRUEorFALSE(wao),
              isTRUEorFALSE(u),
              isTRUEorFALSE(c),
              isTRUEorFALSE(v),
              isSingleNumber(f),
              isTRUEorFALSE(F),
              isTRUEorFALSE(r),
              isTRUEorFALSE(e),
              isTRUEorFALSE(s),
              isTRUEorFALSE(S),
              isTRUEorFALSE(split),
              isTRUEorFALSE(g),
              isTRUEorFALSE(header),
              is.character(names), !anyNA(names), length(names) == length(b),
              isTRUEorFALSE(filenames),
              isTRUEorFALSE(sortout))

    a <- normA(a)
    b <- normB(b)
    R(gr_a <- import(a))
    if (is.character(b) || b[[1L]] == quote(BEDFile))
        R(gr_b <- import(b))
    else {
        R(bv <- b)
        R(grl_b <- List(lapply(bv, import)))
        if (wb || sortout) {
            if (filenames) {
                R(names(grl_b) <- vapply(bv, as.character, character(1L)))
            } else if (!missing(names)) {
                R(names(grl_b) <- names)
            }
            R(gr_b <- stack(grl_b, "b"))
        } else {
            R(gr_b <- unlist(grl_b))
        }
    }

    if (bed) {
        ### FIXME: cannot access gr_a or gr_b!!!
        if (is(gr_a, "GAlignments")) {
            R(gr_a <- grglist(gr_a))
            if (split)
                R(gr_a <- unlist(gr_a, use.names=FALSE))
        }
        if (is(gr_b, "GAlignments")) {
            R(gr_b <- grglist(gr_b))
            if (split)
                R(gr_b <- unlist(gr_b, use.names=FALSE))
        }
    }

    if (wo || wao || loj) {
        wa <- wb <- TRUE
    }
    if (wao) {
        loj <- TRUE
    }

    if (S) {
        gr_b_o <- .R(invert(gr_b))
    } else {
        gr_b_o <- .R(gr_b)
    }

    fracRestriction <- !(missing(f) && missing(F))
    if (!fracRestriction) {
        if (u) {
            return(R(ans <- subsetByOverlaps(gr_a, gr_b_o, ignore.strand=!s)))
        } else if (c) {
            R(ans <- gr_a)
            return(R(mcols(ans)$c <- countOverlaps(gr_a, gr_b_o,
                                                   ignore.strand=!s)))
        } else if (v) {
            if (!s && !S) {
                gr_b_o <- .R(unstrand(gr_b))
            }
            return(R(ans <- gr_a[gr_a %outside% gr_b_o]))
        }
    }
    
    R(hits <- findOverlaps(gr_a, gr_b_o, ignore.strand=!s))

    if (fracRestriction) {
        R(o <- pintersect(gr_a[queryHits(hits)], gr_b[subjectHits(hits)]))
        if (r) {
            F <- f
        }
        if (!missing(f)) {
            keep_f <- .R(o / width(gr_a)[queryHits(hits)] >= f)
        }
        if (!missing(F)) {
            keep_F <- .R(o / width(gr_b)[subjectHits(hits)] >= F)
            if (!missing(f)) {
                keep <- if (e) .R(keep_f | keep_F) else .R(keep_f & keep_F)
            } else {
                keep <- keep_F
            }
        } else {
            keep <- keep_f
        }
        R(hits <- hits[keep])
    }
    
    if (loj) {
        R(seqlevels(gr_b) <- union(".", seqlevels(gr_b)))
        R(ans <- merge(gr_a, gr_b, hits, all.x=loj, y.name="b",
                       NA.VALUE=GRanges(".", IRanges(0L, -1L))))
    } else if (wa && wb) {
        R(ans <- merge(gr_a, gr_b, hits, y.name="b"))
    } else if (wa) {
        R(ans <- gr_a[queryHits(hits)])
    } else {
        R(ans <- pintersect(gr_a[queryHits(hits)], gr_b[subjectHits(hits)]))
        if (wb) {
            R(mcols(ans)$b <- gr_b[subjectHits(hits)])
        }
    }
    
    if (wao) {
        R(seqlevels(gr_a) <- union(".", seqlevels(gr_a)))
    }
    if (wo || wao) {
        R(mcols(ans)$o <- width(pintersect(ans, mcols(ans)$b)))
    }

    if (sortout) {
        R(ans <- sort(ans))
    }
}

BEDTOOLS_INTERSECT_DOC <-
    "Usage:
       bedtools_intersect [options]
     Options:
       -a <FILE>  BAM/BED/GFF/VCF file A. Each feature in A is compared to B in search of overlaps. Use 'stdin' if passing A with a UNIX pipe.
       -b <FILE1>... One or more BAM/BED/GFF/VCF file(s) B. Use 'stdin' if passing B with a UNIX pipe. -b may be followed with multiple databases and/or wildcard (*) character(s).
       -ubam  Write uncompressed BAM output. The default is write compressed BAM output.
       -bed  When using BAM input (-abam), write output as BED. The default is to write output in BAM when using -abam.
       -wa  Write the original entry in A for each overlap.
       -wb  Write the original entry in B for each overlap. Useful for knowing what A overlaps. Restricted by -f and -r.
       -loj  Perform a 'left outer join'. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B.
       -wo  Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r.
       -wao  Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0. Restricted by -f and -r.
       -u  Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B. Restricted by -f and -r.
       -c  For each entry in A, report the number of hits in B while restricting to -f. Reports 0 for A entries that have no overlap with B. Restricted by -f and -r.
       -v  Only report those entries in A that have no overlap in B. Restricted by -f and -r.
       -f  Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
       -F  Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).
       -r  Require that the fraction of overlap be reciprocal for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
       -e  Require that the minimum fraction be satisfied for A _OR_ B. In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of B is covered. Without -e, both fractions would have to be satisfied.
       -s  Force strandedness. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
       -S  Require different strandedness. That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.
       -split  Treat split BAM (i.e., having an 'N' CIGAR operation) or BED12 entries as distinct BED intervals.
       -header  Print the header from the A file prior to results.
       -names <name>... When using multiple databases (-b), provide an alias for each that will appear instead of a fileId when also printing the DB record.
       -filenames  When using multiple databases (-b), show each complete filename instead of a fileId when also printing the DB record.
       -sortout  When using multiple databases (-b), sort the output DB hits for each record."

