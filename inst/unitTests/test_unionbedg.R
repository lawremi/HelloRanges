### =========================================================================
### Tests for bedtools unionbedg command
### -------------------------------------------------------------------------
###
### Based on examples from bedtools (C) 2016 Aaron Quinlan et al.
###

test_unionbedg <- function() {
    setwd(system.file("unitTests", "data", "unionbedg", package="HelloRanges"))

    rowData <-
        GRanges("chr1",
                IRanges(c(901, 1001, 1501, 1701, 1981, 2001, 2051, 2071, 2091),
                        c(1000, 1500, 1600, 1980, 2000, 2050, 2070, 2090, 2100))
                )
    score <- matrix(c(NA, 10, NA, NA, NA, 20, 20, 20, 20, 60, 60, 60,
                      50, 50, 50, NA, NA, NA, NA, NA, NA, NA,
                      80, 80, 80, NA, 20),
                    ncol=3, dimnames=list(NULL, c("1", "2", "3")))
    exp <- SummarizedExperiment(list(score=score), rowData)
    r <- bedtools_unionbedg("-i a.bedGraph,b.bedGraph,c.bedGraph")
    checkEquals(exp, eval(r))

    colnames(exp) <- c("A", "B", "C")
    r <- bedtools_unionbedg("-i a.bedGraph,b.bedGraph,c.bedGraph -names A,B,C")
    checkEquals(exp, eval(r))

    rowData <- GRanges("chr1", IRanges(c(1, 1601, 2101), c(900, 1700, 5000)))
    score <- assay(exp)[rep(NA_integer_, 3),]
    none <- SummarizedExperiment(list(score=score), rowData)
    exp <- sort(rbind(exp, none))
    seqinfo(exp) <- import("test.genome")
    r <- bedtools_unionbedg("-i a.bedGraph,b.bedGraph,c.bedGraph -names A,B,C -empty -g test.genome")
    checkEquals(exp, eval(r))
}
