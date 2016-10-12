### =========================================================================
### Tests for bedtools getfasta command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_getfasta <- function() {
    setwd(system.file("unitTests", "data", "getfasta", package="HelloRanges"))

    fa <- readDNAStringSet("t.fa")
    blocks <- import("blocks.bed")
    exp <- fa[unstrand(blocks)]
    r <- bedtools_getfasta("--fi t.fa -bed blocks.bed")
    checkIdentical(exp, eval(r))

    exp <- fa[unstrand(blocks(blocks))]
    r <- bedtools_getfasta("--fi t.fa -bed blocks.bed -split")
    checkIdentical(exp, eval(r))

    exp <- fa[blocks(blocks)]
    r <- bedtools_getfasta("--fi t.fa -bed blocks.bed -s -split")
    checkIdentical(exp, eval(r))
}
