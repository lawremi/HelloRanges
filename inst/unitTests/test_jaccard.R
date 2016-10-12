### =========================================================================
### Tests for bedtools jaccard command
### -------------------------------------------------------------------------
###
### Based on tests from bedtools (C) 2016 Aaron Quinlan et al.
###

test_jaccard <- function() {
    setwd(system.file("unitTests", "data", "jaccard", package="HelloRanges"))

    exp <- DataFrame(intersection=110L, union=110L, jaccard=1,
                     n_intersections=2L)
    r <- bedtools_jaccard("-a a.bed -b a.bed")
    checkIdentical(exp, eval(r))

    exp <- DataFrame(intersection=5L, union=35L, jaccard=5/35,
                     n_intersections=1L)
    r <- bedtools_jaccard("-a three_blocks_match.bed -b e.bed -split")
    checkIdentical(exp, eval(r))

    exp <- DataFrame(intersection=10L, union=150L, jaccard=10/150,
                     n_intersections=1L)
    r <- bedtools_jaccard("-a a.bam -b three_blocks_match.bam")
    checkIdentical(exp, eval(r))

    exp <- DataFrame(intersection=145L, union=180L, jaccard=145/180,
                     n_intersections=2L)
    r <- bedtools_jaccard("-a aMixedStrands.bed -b bMixedStrands.bed")
    checkIdentical(exp, eval(r))
    
    exp <- DataFrame(intersection=150, union=395, jaccard=150/395,
                     n_intersections=5L)
    r <- bedtools_jaccard("-a aMixedStrands.bed -b bMixedStrands.bed -s")
    checkEquals(exp, eval(r))

    exp <- DataFrame(intersection=70, union=395, jaccard=70/395,
                     n_intersections=2L)
    r <- bedtools_jaccard("-a aMixedStrands.bed -b bMixedStrands.bed -s -f 0.8")
    checkEquals(exp, eval(r))
}
