### =========================================================================
### Tests for bedtools shift command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_shift <- function() {
    setwd(system.file("unitTests", "data", "shift", package="HelloRanges"))

    genome <- import("tiny.genome")
    gr <- import("a.bed", genome=genome)
    exp <- shift(gr, 5L)
    r <- bedtools_shift("-i a.bed -s 5 -g tiny.genome")
    checkIdentical(exp, eval(r))

    r <- bedtools_shift("-i a.bed -p 5 -m 5 -g tiny.genome")
    checkIdentical(exp, eval(r))

    exp <- gr
    exp[strand(exp) == "-"] <- shift(exp[strand(exp) == "-"], 5L)
    r <- bedtools_shift("-i a.bed -p 0 -m 5 -g tiny.genome")
    checkIdentical(exp, eval(r))

    exp <- gr
    exp[strand(exp) == "+"] <- shift(exp[strand(exp) == "+"], 5L)
    r <- bedtools_shift("-i a.bed -p 5 -m 0 -g tiny.genome")
    checkIdentical(exp, eval(r))
}
