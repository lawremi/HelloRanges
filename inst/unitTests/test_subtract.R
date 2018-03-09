### =========================================================================
### Tests for bedtools subtract command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_subtract <- function() {
    setwd(system.file("unitTests", "data", "subtract", package="HelloRanges"))

    a <- import("a.bed")
    b <- import("b.bed")

    hits <- findOverlaps(a, b, ignore.strand = TRUE)
    toSubtract <- reduce(extractList(b, as(hits, "List")),
                         ignore.strand=TRUE)
    exp <- subset(psetdiff(a, toSubtract, ignore.strand=TRUE), width > 0L)
    r <- bedtools_subtract("-a a.bed -b b.bed")
    checkIdentical(exp, eval(r))

    r <- bedtools_subtract("-a a.bed -b b.bed -f 0.1")
    checkIdentical(exp, eval(r))

    exp <- unname(split(granges(a), seq_along(a)))
    r <- bedtools_subtract("-a a.bed -b b.bed -f 0.5")
    checkIdentical(exp, eval(r))

    r <- bedtools_subtract("-a a.bed -b b.bed -s")
    checkIdentical(exp, eval(r))

    exp <- a[2L]
    r <- bedtools_subtract("-a a.bed -b b.bed -A")
    checkIdentical(exp, eval(r))

    r <- bedtools_subtract("-a a.bed -b b.bed -A -f 0.1")
    checkIdentical(exp, eval(r))

    exp <- a
    r <- bedtools_subtract("-a a.bed -b b.bed -A -f 0.5")
    checkIdentical(exp, eval(r))

    c <- import("c.bed")
    exp <- c
    r <- bedtools_subtract("-a c.bed -b d.bed -N -f 0.4")
    checkIdentical(exp, eval(r))

    exp <- c[NULL]
    r <- bedtools_subtract("-a c.bed -b d.bed -N -f 0.39")
    checkIdentical(exp, eval(r))
}
