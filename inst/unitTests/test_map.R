### =========================================================================
### Tests for bedtools map command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_map <- function() {
    setwd(system.file("unitTests", "data", "map", package="HelloRanges"))

    a <- import("ivls.bed")
    b <- import("values.bed")
    hits <- findOverlaps(a, b, ignore.strand=TRUE)
    exp <- a
    mcols(exp) <- aggregate(b, hits, score.sum = sum(score), drop=FALSE)
    
    r <- bedtools_map("-a ivls.bed -b values.bed")
    checkIdentical(exp, eval(r))

    r <- bedtools_map("-a ivls.bed -b values.bed -o sum")
    checkIdentical(exp, eval(r))

    mcols(exp) <- aggregate(b, hits, score.mode = distmode(score), drop=FALSE)
    r <- bedtools_map("-a ivls.bed -b values.bed -o mode")
    checkIdentical(exp, eval(r))

    gff <- import("test.gff2")
    hits <- findOverlaps(a, gff, ignore.strand=TRUE)
    mcols(exp) <- aggregate(gff, hits,
                            seqnames.collapse = unstrsplit(seqnames, ","),
                            drop=FALSE)
    r <- bedtools_map("-a ivls.bed -b test.gff2 -c 1 -o collapse")
    checkIdentical(exp, eval(r))

    a <- import("ivls2.bed")
    b <- import("values5.bed")
    hits <- findOverlaps(a, b, ignore.strand=TRUE)[-1L]
    exp <- a
    mcols(exp) <- aggregate(b, hits,
                            name.collapse = unstrsplit(name, ","),
                            drop=FALSE)
    r <- bedtools_map("-a ivls2.bed -b values5.bed -c 4 -o collapse -f 0.7")
    checkIdentical(exp, eval(r))
}
