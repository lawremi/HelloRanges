### =========================================================================
### Tests for bedtools genomecov command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_genomecov <- function() {
    setwd(system.file("unitTests", "data", "genomecov", package="HelloRanges"))
    
    genome <- import("test.genome")
    y <- import("y.bed", genome=genome)
    cov_gr <- as(coverage(granges(y)), "GRanges")
    exp <- subset(cov_gr, score > 0L)    
    r <- bedtools_genomecov("-i y.bed -bg -g test.genome")
    checkIdentical(exp, eval(r))

    asBam("three_blocks.sam", "three_blocks", overwrite=TRUE)

    three_blocks <- import("three_blocks.bam")
    cov_gr <- as(coverage(granges(three_blocks)), "GRanges")
    exp <- subset(cov_gr, score > 0L)
    r <- bedtools_genomecov("-i three_blocks.bam -bg")
    checkIdentical(exp, eval(r))

    cov_gr <- as(coverage(three_blocks), "GRanges")
    exp <- subset(cov_gr, score > 0L)
    r <- bedtools_genomecov("-i three_blocks.bam -bg -split")
    checkIdentical(exp, eval(r))

    exp <- cov_gr
    r <- bedtools_genomecov("-i three_blocks.bam -bga -split")
    checkIdentical(exp, eval(r))

    exp <- GPos(coverage(three_blocks))
    r <- bedtools_genomecov("-i three_blocks.bam -dz -split")
    checkIdentical(exp, eval(r))

    asBam("sam-w-del.sam", "sam-w-del", overwrite=TRUE)
    sam_w_del <- import("sam-w-del.bam")

    cov_gr <- as(coverage(sam_w_del, drop.D.ranges=TRUE), "GRanges")
    exp <- subset(cov_gr, score > 0L)
    r <- bedtools_genomecov("-i sam-w-del.bam -bg -split")
    checkIdentical(exp, eval(r))

    exp <- DataFrame(seqnames=Rle(c("1", "2", "3", "genome"),
                                  c(3L, 1L, 1L, 3L)),
                     coverage = factor(c(0:2, 0, 0, 0:2)),
                     count = c(93, 4, 3, rep(100, 2), 293, 4, 3),
                     len = c(rep(100, 5), rep(300, 3)))
    exp <- within(exp, fraction <- count / len)
    r <- bedtools_genomecov("-i y.bam")
    checkIdentical(exp, eval(r))

    asBam("pair-chip.sam", "pair-chip", overwrite=TRUE)
    pair_chip <- import("pair-chip.bam", paired=TRUE)
    cov_gr <- as(coverage(granges(pair_chip)), "GRanges")
    exp <- subset(cov_gr, score > 0L)
    r <- bedtools_genomecov("-i pair-chip.bam -bg -pc")
    checkIdentical(exp, eval(r))

    asBam("chip.sam", "chip", overwrite=TRUE)
    chip <- import("chip.bam")
    cov_gr <- as(coverage(resize(granges(chip), 100L)), "GRanges")
    exp <- subset(cov_gr, score > 0L)
    r <- bedtools_genomecov("-i chip.bam -bg -fs 100")
    checkIdentical(exp, eval(r))
}
