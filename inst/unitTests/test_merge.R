### =========================================================================
### Tests for bedtools merge command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_merge <- function() {
    setwd("data/merge")

    a <- import("a.bed")

    exp <- reduce(a)
    r <- bedtools_merge("-i a.bed")
    checkIdentical(exp, eval(r))

    exp <- reduce(a, with.revmap=TRUE)
    mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                            seqnames.count = lengths(seqnames))
    r <- bedtools_merge("-i a.bed -c 1 -o count")
    checkIdentical(exp, eval(r))

    mcols(a)$name <- paste0("a", 1:4)
    exp <- reduce(a, with.revmap=TRUE)
    mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                            name.collapse = unstrsplit(name, ","))
    r <- bedtools_merge("-i a.names.bed -c 4  -o collapse")
    checkIdentical(exp, eval(r))

    exp <- reduce(a, with.revmap=TRUE)
    mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                            name.collapse = unstrsplit(name, "|"))
    r <- bedtools_merge("-i a.names.bed -delim \"|\" -c 4  -o collapse")
    checkIdentical(exp, eval(r))

    a <- import("a.full.bed")
    exp <- reduce(a, with.revmap=TRUE, ignore.strand=TRUE)
    mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                            name.collapse = unstrsplit(name, ","),
                            score.sum = sum(score))
    r <- bedtools_merge("-i a.full.bed -c 4,5  -o collapse,sum")
    checkIdentical(exp, eval(r))

    exp <- reduce(a, with.revmap=TRUE, ignore.strand=TRUE)
    mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                            score.count = lengths(score),
                            score.sum = sum(score))
    r <- bedtools_merge("-i a.full.bed -c 5  -o count,sum")
    checkIdentical(exp, eval(r))

### FIXME: we should not need to tabix the VCF, but VcfFile (and
### thus import) requires it!
    p <- bgzip("testA.vcf", overwrite=TRUE)
    indexTabix(p, format="vcf", overwrite=TRUE)
    vcf <- readVcf(p, genome=Seqinfo())
    exp <- reduce(granges(vcf))
    r <- bedtools_merge("-i testA.vcf.bgz")
    checkIdentical(exp, eval(r))

    a <- import("a.full.bed")
    exp <- reduce(subset(a, strand=="+"))
    r <- bedtools_merge("-i a.full.bed -S +")
    checkIdentical(exp, eval(r))

    bam <- import("fullFields.bam", use.names=TRUE)
    exp <- reduce(granges(bam), with.revmap=TRUE, ignore.strand=TRUE)
    mcols(exp) <- aggregate(bam, mcols(exp)$revmap,
                            names.collapse = unstrsplit(names, ","))
    r <- bedtools_merge("-i fullFields.bam -c 1 -o collapse")
    checkIdentical(exp, eval(r))

    exp <- reduce(granges(bam), with.revmap=TRUE, ignore.strand=TRUE)
    mcols(exp) <- aggregate(bam, mcols(exp)$revmap,
                            seqnames.collapse = unstrsplit(seqnames, ","))
    r <- bedtools_merge("-i fullFields.bam -c 3 -o collapse")
    checkIdentical(exp, eval(r))

    exp <- reduce(granges(bam), with.revmap=TRUE, ignore.strand=TRUE)
    mcols(exp) <- aggregate(bam, mcols(exp)$revmap, start.mean = mean(start))
    r <- bedtools_merge("-i fullFields.bam -c 4 -o mean")
    checkIdentical(exp, eval(r))

    bam <- import("fullFields.bam", param = ScanBamParam(what="mapq"))
    exp <- reduce(granges(bam), with.revmap=TRUE, ignore.strand=TRUE)
    mcols(exp) <- aggregate(bam, mcols(exp)$revmap, mapq.mean = mean(mapq))
    r <- bedtools_merge("-i fullFields.bam -c 5 -o mean")
    checkIdentical(exp, eval(r))  
}
