### =========================================================================
### Tests for bedtools closest command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_closest <- function() {
    setwd("data/closest")

    a <- import("a.bed")
    b <- import("b.bed")

    exp <- Pairs(a, b)
    mcols(exp)$distance <- 1L
    r <- bedtools_closest("-a a.bed -b b.bed -d")
    checkIdentical(exp, eval(r))
}
