### =========================================================================
### Tests for bedtools intersect command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_basic_self_intersection <- function() {
    exp <- GRanges("chr1", IRanges(c(10, 100), c(20, 200)), name=c("a1", "a2"),
                   score=c(1, 2), strand=c("+", "-"))
    r <- bedtools_intersect("a.bed", "b.bed")
    checkIdentical(exp, eval(r))
}
