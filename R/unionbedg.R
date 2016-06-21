### =========================================================================
### bedtools unionbedg command
### -------------------------------------------------------------------------
###

bedtools_unionbedg <- function(cmd = "--help") {
    do_R_call(R_bedtools_unionbedg, BEDTOOLS_UNIONBEDG_DOC, cmd)
}

R_bedtools_unionbedg <- function(i, header=FALSE, names=NULL,
                                 g=NA_character_, empty=FALSE)
{
    .R_bedtools_disjoin(i, header, names, g, empty, use.score=TRUE)
}

BEDTOOLS_UNIONBEDG_DOC <-
    "Usage:
       bedtools_unionbedg [options]
     Options:
           -i <FILE,...>  bedGraph files.
      -header  Print a header line. (chrom/start/end + names of each file).
      --names <name,...>  A list of names to describe each file in -i.
               These names will be printed in the header line.
           -g <id>  Use genome file to calculate empty regions.
      --empty  Report empty regions (i.e., start/end intervals w/o
               values in all files)."

do_bedtools_unionbedg <- make_do(R_bedtools_unionbedg)
