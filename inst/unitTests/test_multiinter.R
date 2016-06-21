### =========================================================================
### Tests for bedtools multiinter command
### -------------------------------------------------------------------------
###
### Based on examples from bedtools (C) 2016 Aaron Quinlan et al.
###

test_multiinter <- function() {
    setwd(system.file("unitTests", "data", "multiinter", package="HelloRanges"))
    
    exp <- GRanges("chr1",
                   IRanges(c(7, 9, 13, 16, 21, 23, 31, 33),
                           c(8, 12, 15, 20, 22, 30, 32, 34)),
                   revmap=IntegerList(1, c(1,4), c(1,3,4), c(1,3), 3,
                                      c(2,3), 3, 5),
                   i=FactorList(1, c(1,3), 1:3, 1:2, 2, 1:2, 2, 3))
    r <- bedtools_multiinter("-i a.bed,b.bed,c.bed")
    checkIdentical(exp, eval(r))

    mcols(exp)$i <- extractList(factor(c("A", "B", "C")),
                                IntegerList(mcols(exp)$i))
    r <- bedtools_multiinter("-i a.bed,b.bed,c.bed -names A,B,C")    
    checkIdentical(exp, eval(r))

    none <- GRanges("chr1", IRanges(c(1, 35), c(6, 5000)))
    mcols(none)$revmap <- IntegerList(integer())
    mcols(none)$i <- FactorList(factor())
    exp <- sort(c(exp, none))
    seqinfo(exp) <- import("test.genome")
    r <- bedtools_multiinter("-i a.bed,b.bed,c.bed -names A,B,C -empty -g test.genome")
    checkIdentical(exp, eval(r))
}
