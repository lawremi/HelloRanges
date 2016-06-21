### =========================================================================
### utilities for defining the package
### -------------------------------------------------------------------------
###

make_do <- function(fun) {
    args <- formals(fun)
    args[] <- lapply(names(args), as.name)
    as.function(c(formals(fun),
                  call("eval", as.call(c(substitute(fun), args)))))
}

globalVariables(c(".gr_a_o", ".gr_b_o", "BEDFile", "ScanBamParam",
                  "SummarizedExperiment", "alphabetFrequency", "asBED",
                  "count", "genome", "ignore.strand", "keep", "len",
                  "letterFrequency", "olap", "pairs", "readDNAStringSet",
                  "seqlengths", "vcountPattern"))

