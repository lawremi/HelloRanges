### =========================================================================
### Tests for bedtools groupby command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_groupby <- function() {
    setwd(system.file("unitTests", "data", "groupby", package="HelloRanges"))

    a <- import("values3.header.bed")
    exp <- aggregate(unstrand(a), score.sum = sum(score))
    r <- bedtools_groupby("-i values3.header.bed -c 5")
    checkIdentical(exp, eval(r))

    indexTabix(bgzip("a_vcfSVtest.vcf", overwrite=TRUE), "vcf")
    a <- import("a_vcfSVtest.vcf.bgz")
    exp <- aggregate(granges(a), ~seqnames + REF, QUAL.mean = mean(QUAL))
    r <- bedtools_groupby("-i a_vcfSVtest.vcf.bgz -g 1,4 -c 6 -o mean")
    checkIdentical(exp, eval(r))
}
