### =========================================================================
### Tests for bedtools makewindows command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_makewindows <- function() {
    setwd("data/makewindows")

    input <- import("input.bed")

    exp <- tile(input, w=5000L)
    r <- bedtools_makewindows("-b input.bed -w 5000")
    checkIdentical(exp, eval(r))

    exp <- slidingWindows(input, w=5000L, s=2000L)
    r <- bedtools_makewindows("-b input.bed -w 5000 -s 2000")
    checkIdentical(exp, eval(r))

    exp <- tile(input, 3L)
    r <- bedtools_makewindows("-b input.bed -n 3")
    checkIdentical(exp, eval(r))

    genome <- import("test.genome")
    exp <- tile(as(genome, "GRanges"), 3L)
    r <- bedtools_makewindows("-g test.genome -n 3")
    checkIdentical(exp, eval(r))
}
