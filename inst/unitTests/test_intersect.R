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

    r <- bedtools_intersect("-a a.bed -b a.bed -v")
    checkIdentical(exp[NULL], eval(r))

    cexp <- exp
    mcols(cexp)$c <- c(0L, 2L)
    r <- bedtools_intersect("-a a.bed -b b.bed -c")
    checkIdentical(cexp, eval(r))

    mcols(cexp)$c <- c(0L, 1L)
    r <- bedtools_intersect("-a a.bed -b b.bed -c -s")
    checkIdentical(cexp, eval(r))

    mcols(cexp)$c <- c(0L, 0L)
    r <- bedtools_intersect("-a a.bed -b b.bed -c -s -f 0.1")
    checkIdentical(cexp, eval(r))

    exp <- GRanges("chr1", IRanges(c(101, 101), c(101, 110)),
                   name=c("a2", "a2"),
                   score=c(2, 2), strand=c("-", "-"))
    r <- bedtools_intersect("-a a.bed -b b.bed")
    checkIdentical(exp, fixup(eval(r)))

    exp_a <- import("a.bed")
    exp_b <- import("b.bed")
    
    r <- bedtools_intersect("-a a.bed -b b.bed -wa")
    checkIdentical(exp_a[c(2, 2)], eval(r))

    exp_a_b <- Pairs(exp_a[c(2, 2)], exp_b[2:3])

    r <- bedtools_intersect("-a a.bed -b b.bed -wa -wb")
    checkIdentical(exp_a_b, eval(r))

    exp_o <- exp_a_b
    mcols(exp_o)$o <- c(1L, 10L)
    r <- bedtools_intersect("-a a.bed -b b.bed -wa -wb -wo")
    checkIdentical(exp_o, eval(r))

    suppressWarnings({
        first <- exp_a[c(1, 2, 2)]
        seqlevels(first) <- c(".", seqlevels(first))
        exp_loj <- Pairs(first, c(NAGRanges(exp_b), exp_b[2:3]))
    })
    mcols(exp_loj)$o <- c(0L, 1L, 10L)

    r <- bedtools_intersect("-a a.bed -b b.bed -wa -wb -wao")
    checkIdentical(exp_loj, eval(r))

}
