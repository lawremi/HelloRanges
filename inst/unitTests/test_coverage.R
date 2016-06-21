### =========================================================================
### Tests for bedtools coverage command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_coverage <- function() {
    setwd(system.file("unitTests", "data", "coverage", package="HelloRanges"))
    
    a <- import("a.bed")
    b <- import("b.bed")
    genome <- import("test.genome")
    
    exp <- a
    mcols(exp) <- cbind(mcols(exp),
                        DataFrame(count=c(2L, 5L, 4L, 6L, 7L, 6L),
                                  covered=c(50L, 50L, 38L, 50L, 50L, 50L)))
    mcols(exp)$fraction <- with(exp, covered/width)
    r <- bedtools_coverage("-a a.bed -b b.bed")
    checkIdentical(exp, eval(r))

    mcols(exp) <- mcols(exp)[1:3]
    r <- bedtools_coverage("-a a.bed -b b.bed -counts")
    checkIdentical(exp, eval(r))

    coverage <- factor(c(2, 2, 3, 0, 1:4, 3:6, 1:6, 0:6))
    count <- c(50, 40, 10, 12, 20, 11, 4, 3, 46, 4, 22, 28, 16, 2, 6, 4, 13, 9,
               12, 36, 103, 66, 11, 35, 37)
    p <- PartitioningByWidth(c(1, 2, 5, 2, 2, 6, 7))
    len <- rep(rep(c(50L, 300L), c(6, 1)), width(p))
    covhist <- DataFrame(coverage, count, len, fraction=count/len)
    covhistList <- relist(covhist, p)
    names(covhistList) <- c(1:6, "all")
    exp <- a
    seqinfo(exp) <- import("test.genome")
    mcols(exp)$coverage <- head(covhistList, -1L)
    metadata(exp)$coverage <- covhistList$all
    r <- bedtools_coverage("-a a.bed -b b.bed -hist -g test.genome")
    checkIdentical(exp, eval(r))

    exp <- a
    seqinfo(exp) <- genome
    seqinfo(b) <- genome
    mcols(exp)$coverage <- unname(coverage(b)[exp])
    r <- bedtools_coverage("-a a.bed -b b.bed -d -g test.genome")
    checkIdentical(exp, eval(r))

    cov_pos <- coverage(subset(b, strand=="+"))[subset(exp, strand=="+")]
    cov_neg <- coverage(subset(b, strand=="-"))[subset(exp, strand=="-")]
    cov_star <- coverage(b)[subset(exp, strand=="*")]
    mcols(exp)$coverage <- unname(unsplit(List(cov_pos, cov_neg, cov_star),
                                          strand(exp)))
    r <- bedtools_coverage("-a a.bed -b b.bed -d -s -g test.genome")
    checkIdentical(exp, eval(r))
    
    exp <- a
    mcols(exp) <- cbind(mcols(exp),
                        DataFrame(count=c(2L, 1L, 0L, 4L, 3L, 4L),
                                  covered=c(50L, 23L, 0L, 50L, 50L, 34L)))
    mcols(exp)$fraction <- with(exp, covered/width)
    r <- bedtools_coverage("-a a.bed -b b.bed -s")
    checkIdentical(exp, eval(r))

    exp <- a
    mcols(exp) <- cbind(mcols(exp),
                        DataFrame(count=c(0L, 4L, 4L, 2L, 4L, 2L),
                                  covered=c(0L, 50L, 38L, 50L, 50L, 50L)))
    mcols(exp)$fraction <- with(exp, covered/width)
    r <- bedtools_coverage("-a a.bed -b b.bed -S")
    checkIdentical(exp, eval(r))
}
