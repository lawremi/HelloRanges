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
    mcols(cexp)$overlap_count <- c(0L, 2L)
    r <- bedtools_intersect("-a a.bed -b b.bed -c")
    checkIdentical(cexp, eval(r))

    mcols(cexp)$overlap_count <- c(0L, 1L)
    r <- bedtools_intersect("-a a.bed -b b.bed -c -s")
    checkIdentical(cexp, eval(r))

    mcols(cexp)$overlap_count <- c(0L, 0L)
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
    mcols(exp_o)$overlap_width <- c(1L, 10L)
    r <- bedtools_intersect("-a a.bed -b b.bed -wa -wb -wo")
    checkIdentical(exp_o, eval(r))

    suppressWarnings({
        first <- exp_a[c(1, 2, 2)]
        seqlevels(first) <- c(".", seqlevels(first))
        exp_loj <- Pairs(first, c(NAGRanges(exp_b), exp_b[2:3]))
    })
    mcols(exp_loj)$overlap_width <- c(0L, 1L, 10L)

    r <- bedtools_intersect("-a a.bed -b b.bed -wa -wb -wao")
    checkIdentical(exp_loj, eval(r))

    r <- bedtools_intersect("-a a.bed -b b.bed -wa -wb -wo -s")
    checkIdentical(exp_o[1L], eval(r))

    r <- bedtools_intersect("-a a.bed -b b.bed -wa -wb -wao -s")
    checkIdentical(exp_loj[1:2], eval(r))

    p <- pipe("cat a.bed | Rscript -e 'library(HelloRanges); export(eval(bedtools_intersect(\"-a stdin -b b.bed\")), stdout(), format=\"bed\")'", "r")
    checkIdentical(exp, import(p, format="bed"))
    close(p)

    p <- pipe("cat b.bed | Rscript -e 'library(HelloRanges); export(eval(bedtools_intersect(\"-a a.bed -b stdin\")), stdout(), format=\"bed\")'", "r")
    checkIdentical(exp, import(p, format="bed"))
    close(p)

    one_block <- Rsamtools::asBam("one_block.sam", "one_block", overwrite=TRUE)
    two_blocks <- Rsamtools::asBam("two_blocks.sam", "two_blocks",
                                   overwrite=TRUE)
    three_blocks <- Rsamtools::asBam("three_blocks.sam", "three_blocks",
                                     overwrite=TRUE)

    three_blocks_exp <- GenomicAlignments::readGAlignments(three_blocks)
    r <- bedtools_intersect("-a three_blocks.bam -b three_blocks_nomatch.bed")
    checkIdentical(three_blocks_exp, eval(r))

    r <- bedtools_intersect("-a three_blocks.bam -b three_blocks_nomatch.bed -split")
    checkIdentical(three_blocks_exp[NULL], eval(r))
    
    r <- bedtools_intersect("-a three_blocks.bam -b three_blocks_match.bed -split")
    checkIdentical(three_blocks_exp, eval(r))

    r <- bedtools_intersect("-a three_blocks.bam -b three_blocks_match.bed -split -s")
    checkIdentical(three_blocks_exp[NULL], eval(r))

    r <- bedtools_intersect("-a three_blocks.bam -b three_blocks_match_1bp.bed -split -f 0.1")
    checkIdentical(three_blocks_exp[NULL], eval(r))

    one_block_c_exp <- GenomicAlignments::readGAlignments(one_block)
    mcols(one_block_c_exp)$overlap_count <- 1L
    r <- bedtools_intersect("-a one_block.bam -b c.bed -c")
    checkIdentical(one_block_c_exp, eval(r))

    one_block_exp <- GenomicAlignments::readGAlignments(one_block)
    exp_c <- import("c.bed")
    bam_wo_exp <- Pairs(one_block_exp, exp_c)
    mcols(bam_wo_exp)$overlap_width <- 30L
    r <- bedtools_intersect("-a one_block.bam -b c.bed -wo")
    checkIdentical(bam_wo_exp, eval(r))

    seqlevels(first(bam_wo_exp)) <- c(".", seqlevels(first(bam_wo_exp)))
    seqlevels(second(bam_wo_exp)) <- c(".", seqlevels(second(bam_wo_exp)))
    r <- bedtools_intersect("-a one_block.bam -b c.bed -wao")
    checkIdentical(bam_wo_exp, eval(r))
}
