### =========================================================================
### Tests for bedtools slop command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_slop <- function() {
    setwd("data/slop")

    genome <- import("tiny.genome")
    a <- import("a.bed", genome=genome)
    
    exp <- a + 5L
    r <- bedtools_slop("-i a.bed -b 5 -g tiny.genome")
    checkIdentical(exp, eval(r))

    r <- bedtools_slop("-i a.bed -l 5 -r 5 -g tiny.genome")
    checkIdentical(exp, eval(r))

    exp <- resize(a, 5 + width(a), fix="end", ignore.strand=TRUE)
    r <- bedtools_slop("-i a.bed -l 5 -r 0 -g tiny.genome")
    checkIdentical(exp, eval(r))

    exp <- resize(a, 5 + width(a), fix="start", ignore.strand=TRUE)
    r <- bedtools_slop("-i a.bed -l 0 -r 5 -g tiny.genome")
    checkIdentical(exp, eval(r))

    exp <- resize(a, 5 + width(a), fix="start", ignore.strand=FALSE)
    r <- bedtools_slop("-i a.bed -r 5 -s -g tiny.genome")
    checkIdentical(exp, eval(r))
}
