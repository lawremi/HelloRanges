### =========================================================================
### Tests for bedtools flank command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_flank <- function() {
    setwd(system.file("unitTests", "data", "flank", package="HelloRanges"))

    genome <- import("tiny.genome")
    gr_a <- import("a.bed", genome=genome)
    left <- flank(gr_a, 5L, ignore.strand = TRUE)
    right <- flank(gr_a, 5L, start = FALSE, ignore.strand = TRUE)
    exp <- zipup(Pairs(left, right))
    r <- bedtools_flank("-i a.bed -b 5 -g tiny.genome")
    checkIdentical(exp, eval(r))

    r <- bedtools_flank("-i a.bed -l 5 -r 5 -g tiny.genome")
    checkIdentical(exp, eval(r))

    exp <- left
    r <- bedtools_flank("-i a.bed -l 5 -r 0 -g tiny.genome")
    checkIdentical(exp, eval(r))

    exp <- right
    r <- bedtools_flank("-i a.bed -l 0 -r 5 -g tiny.genome")
    checkEquals(exp, eval(r))

    exp <- flank(gr_a, 5L, ignore.strand = FALSE)
    r <- bedtools_flank("-i a.bed -l 5 -r 0 -s -g tiny.genome")
    checkIdentical(exp, eval(r))

    r <- bedtools_flank("-i a.bed -l 0.5 -r 0 -pct -g tiny.genome")
}
