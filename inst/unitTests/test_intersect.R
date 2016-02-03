### =========================================================================
### Tests for bedtools intersect command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

library(HelloRanges)
library(RUnit)

test_basic_self_intersection <- function() {
    fixup <- function(x) {
        mcols(x)$hit <- NULL # 'hit' column added by pintersect()
        x
    }
    setwd("data/intersect")
    exp <- GRanges("chr1", IRanges(c(11, 101), c(20, 200)), name=c("a1", "a2"),
                   score=c(1, 2), strand=c("+", "-"))
    r <- bedtools_intersect("-a a.bed -b a.bed")
    checkIdentical(exp, fixup(eval(r)))
}
