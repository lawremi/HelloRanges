### =========================================================================
### bedtools multiinter command
### -------------------------------------------------------------------------
###

bedtools_multiinter <- function(cmd = "--help") {
    do_R_call(R_bedtools_multiinter, BEDTOOLS_MULTIINTER_DOC, cmd)
}

.R_bedtools_disjoin <- function(i, header=FALSE, names=NULL,
                                g=NA_character_, empty=FALSE, use.score=FALSE)
{
    stopifnot((is.character(i) && !anyNA(i) && length(i) >= 1L) || hasRanges(i),
              isTRUEorFALSE(header),
              isGenome(g),
              isTRUEorFALSE(empty),
              isTRUEorFALSE(use.score))

    if (!is.null(names)) {
        stopifnot(is.character(names), !anyNA(names),
                  length(names) == length(i))
    }

    importGenome(g)

    i <- normB(i)
    if (use.score) {
        .gr_i <- importB(i, names)
    } else {
        .gr_i <- importB(i, names, beforeStack=reduce)
    }
    .gr_i_o <- prepOverlapRanges(i)
    
    R(dj <- disjoin(.gr_i, ignore.strand=TRUE, with.revmap=TRUE))
    rm(i)    
    if (use.score) {
        R(mcols(.gr_i)$i <- decode(mcols(.gr_i)$i))
        R(dfl <- extractList(mcols(.gr_i), mcols(dj)$revmap))
        R(assay <- as.matrix(dfl[,"score"], col.names=dfl[,"i"]))
        R(rowData <- granges(dj))
    } else {
        R(mcols(dj)$i <- extractList(decode(mcols(.gr_i)$i), mcols(dj)$revmap))
    }

    if (empty) {
        R(none <- setdiff(as(seqinfo(dj), "GRanges"), dj))
        if (use.score) {
            R(assay <- rbind(assay,
                             matrix(nrow=length(none), ncol=ncol(assay))))
            R(rowData <- c(rowData, none))
        } else {
            R(mcols(none)$revmap <- IntegerList(integer()))
            R(mcols(none)$i <- FactorList(factor()))
            R(dj <- sort(c(dj, none)))
        }
    }

    if (use.score) {
        R(ans <- SummarizedExperiment(list(score=assay), rowData))
        if (empty) {
            R(ans <- sort(ans))
        }
    } else {
        R(ans <- dj)
    }
    
    R(ans)
}

R_bedtools_multiinter <- function(i, header=FALSE, names=NULL,
                                  g=NA_character_, empty=FALSE)
{
    .R_bedtools_disjoin(i, header, names, g, empty)
}

BEDTOOLS_MULTIINTER_DOC <-
    "Usage:
       bedtools_multiinter [options]
     Options:
           -i <FILE,...>  BAM/BED/GFF/VCF files.
      -header  Print a header line. (chrom/start/end + names of each file).
      --names <name,...>  A list of names to describe each file in -i.
               These names will be printed in the header line.
           -g <id>  Use genome file to calculate empty regions.
      --empty  Report empty regions (i.e., start/end intervals w/o
               values in all files)."

do_bedtools_multiinter <- make_do(R_bedtools_multiinter)

