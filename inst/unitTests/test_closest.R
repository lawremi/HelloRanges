### =========================================================================
### Tests for bedtools closest command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_closest <- function() {
    setwd(system.file("unitTests", "data", "closest", package="HelloRanges"))

    addDotSeq <- function(x) {
        seqlevels(x) <- c(".", seqlevels(x))
        x
    }
    
    a <- addDotSeq(import("a.bed"))
    b <- addDotSeq(import("b.bed"))

    exp <- Pairs(a, b)
    mcols(exp)$distance <- 0L
    r <- bedtools_closest("-a a.bed -b b.bed -d")
    checkIdentical(exp, eval(r))

    exp <- Pairs(b, a)
    mcols(exp)$distance <- 0L
    r <- bedtools_closest("-a b.bed -b a.bed -d")
    checkIdentical(exp, eval(r))

    a <- addDotSeq(import("strand-test-a.bed"))
    b <- addDotSeq(import("strand-test-b.bed"))
    exp <- pair(a, b, nearest(a, b, select="all"), all.x=TRUE)
    r <- bedtools_closest("-a strand-test-a.bed -b strand-test-b.bed -s")
    checkIdentical(exp, eval(r))

    exp <- pair(a, b, nearest(a, b, select="all", ignore.strand=TRUE),
                all.x=TRUE)
    r <- bedtools_closest("-a strand-test-a.bed -b strand-test-b.bed -S")
    checkIdentical(exp, eval(r))

    a <- addDotSeq(import("close-a.bed"))
    b <- addDotSeq(import("close-b.bed"))
    exp <- pair(a, b, breakTies(nearest(a, b, select="all"), "first"),
                all.x=TRUE)
    r <- bedtools_closest("-a close-a.bed -b close-b.bed -t first")
    checkIdentical(exp, eval(r))

    exp <- pair(a, b, breakTies(nearest(a, b, select="all"), "last"),
                all.x=TRUE)
    r <- bedtools_closest("-a close-a.bed -b close-b.bed -t last")
    checkIdentical(exp, eval(r))

    a <- import("mq1.bed")
    mdb1 <- import("mdb1.bed")
    mdb2 <- import("mdb2.bed")
    mdb3 <- import("mdb3.bed")
    b <- mstack(mdb1, mdb2, mdb3, .index.var="b")
    exp <- unlist(List(lapply(split(b, ~ b), function(bi) {
        pair(a, bi, nearest(a, bi, select="all", ignore.strand=TRUE),
             all.x=TRUE)
    })), use.names=FALSE)
    r <- bedtools_closest("-a mq1.bed -b mdb1.bed,mdb2.bed,mdb3.bed")
    checkIdentical(exp, eval(r))

    mcols(second(exp))$b <- Rle(factor(c("a", "b", "c")))
    r <- bedtools_closest("-a mq1.bed -b mdb1.bed,mdb2.bed,mdb3.bed -names a,b,c")
    checkIdentical(exp, eval(r))

    mcols(second(exp))$b <- Rle(factor(c("mdb1.bed", "mdb2.bed", "mdb3.bed")))
    r <- bedtools_closest("-a mq1.bed -b mdb1.bed,mdb2.bed,mdb3.bed -filenames")
    checkIdentical(exp, eval(r))

    b <- mstack(mdb1, mdb2, mdb3, .index.var="b")
    exp <- pair(a, b, nearest(a, b, select="all", ignore.strand=TRUE),
                all.x=TRUE)
    r <- bedtools_closest("-a mq1.bed -b mdb1.bed,mdb2.bed,mdb3.bed -mdb all")
    checkIdentical(exp, eval(r))

    exp <- pair(a, mdb1, nearest(a, mdb1, select="all"), all.x=TRUE)
    r <- bedtools_closest("-a mq1.bed -b mdb1.bed -s")
    checkIdentical(exp, eval(r))

    exp <- pair(a, mdb1, nearest(a, invertStrand(mdb1), select="all"),
                all.x=TRUE)
    r <- bedtools_closest("-a mq1.bed -b mdb1.bed -S")
    checkIdentical(exp, eval(r))


    a <- import("d.bed")
    b <- import("d_iu.bed")
    hits <- precede(a, b, ignore.strand = TRUE, select = "all")
    exp <- pair(a, b, hits, all.x=TRUE)
    mcols(exp)$distance <- distance(exp)
    r <- bedtools_closest("-a d.bed -b d_iu.bed -D ref -iu")
    checkIdentical(exp, eval(r))
    
    b <- import("d_id.bed")
    hits <- follow(a, b, ignore.strand = TRUE, select = "all")
    exp <- pair(a, b, hits, all.x=TRUE)
    mcols(exp)$distance <- distance(exp)
    r <- bedtools_closest("-a d.bed -b d_id.bed -D ref -id")
    checkIdentical(exp, eval(r))
}
