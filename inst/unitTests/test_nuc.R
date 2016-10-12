### =========================================================================
### Tests for bedtools nuc command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_nuc <- function() {
    setwd(system.file("unitTests", "data", "nuc", package="HelloRanges"))

    seq <- readDNAStringSet("test.fasta")
    names(seq) <- "chr1"
    gr <- import("a.bed")
    exp <- seq[gr]
    GC <- as.vector(letterFrequency(exp, "GC", as.prob = TRUE))
    AT <- 1 - GC
    baseCounts <- alphabetFrequency(exp, baseOnly = TRUE)
    mcols(exp) <- DataFrame(AT, GC, baseCounts)

    r <- bedtools_nuc("--fi test.fasta -bed a.bed")
    checkEquals(exp, eval(r))

    mcols(exp)$npattern <- c(0L, 2L)
    r <- bedtools_nuc("--fi test.fasta -bed a.bed -pattern ATA")
    checkEquals(exp, eval(r))
}
