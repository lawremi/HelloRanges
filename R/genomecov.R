### =========================================================================
### bedtools genomecov command
### -------------------------------------------------------------------------
###

bedtools_genomecov <- function(cmd = "--help") {
    do_R_call(R_bedtools_genomecov, BEDTOOLS_GENOMECOV_DOC, cmd)
}

R_bedtools_genomecov <- function(i, g=NA_character_,
                                 d=FALSE, dz=FALSE, bg=FALSE, bga=FALSE,
                                 split=FALSE, strand = c("any", "+", "-"),
                                 `5`=FALSE, `3`=FALSE, max=NULL, scale=1.0,
                                 pc=FALSE, fs=NULL)
{
    strand <- match.arg(strand)
    max <- if (!is.null(max)) as.integer(max)
    fs <- if (!is.null(fs)) as.integer(fs)
    
    stopifnot(isSingleString(i) || hasRanges(i),
              isSingleStringOrNA(g),
              isTRUEorFALSE(d),
              isTRUEorFALSE(dz),
              isTRUEorFALSE(bg),
              isTRUEorFALSE(bga), (d + dz + bg + bga) <= 1L,
              isTRUEorFALSE(split),
              isTRUEorFALSE(`5`),
              isTRUEorFALSE(`3`), !(`5` && `3`),
              is.null(max) || (isSingleNumber(max) && max >= 0L),
              isSingleNumber(scale),
              isTRUEorFALSE(pc), !(pc && !isBam(i)),
              is.null(fs) || (isSingleNumber(fs) && fs > 0L && isBam(i)))

    importGenome(g)

    i <- normA(i)
    if (pc) {
        .gr_i <- importA(i, paired=TRUE)
    } else {
        .gr_i <- importA(i)
    }
    .gr_i_o <- prepOverlapRanges(i, split)
    if (split && isBam(i)) {
        .gr_i_o$drop.D.ranges <- TRUE
    }

    if (strand != "any") {
        .strand <- strand
        rm(strand)
        .gr_i_o <- .R(subset(.gr_i_o, strand == .strand))
    }

    if (!is.null(fs)) {
        .gr_i_o <- .R(resize(.gr_i_o, fs))
    }

    if (`5`) {
        .gr_i_o <- .R(start(.gr_i_o))
    } else if (`3`) {
        .gr_i_o <- .R(end(.gr_i_o))
    }

    R(cov <- coverage(.gr_i_o))

    if (scale != 1.0) {
        R(cov <- cov * scale)
    }
    
    hist <- !(d || dz || bg || bga)
    if (hist) {
        if (!is.null(max)) {
            R(cov[cov > max] <- max)
        }
        R(tablist <- List(lapply(cov, table)))
        R(mcols(tablist)$len <- lengths(cov, use.names=FALSE))
        R(covhist <- stack(tablist, "seqnames", "count", "coverage"))
        R(margin <- aggregate(covhist, ~ coverage,
                              count = sum(NumericList(count)))[-1L])
        R(margin <- DataFrame(seqnames=Rle("genome"), margin,
                              len=sum(as.numeric(lengths(cov)))))
        R(covhist <- rbind(covhist, margin))
        R(ans <- within(covhist, fraction <- count / len))
    } else if (bg || bga) {
        R(ans <- GRanges(cov))
    } else if (d || dz) {
        R(ans <- GPos(cov))
    }

    if (bg || d) {
        R(ans <- subset(ans, score > 0))
    }

    R(ans)
}

BEDTOOLS_GENOMECOV_DOC <-
    "Usage:
       bedtools_genomecov [options]
     Options:
       -i <FILE>  BAM/BED/GFF/VCF file A. Each feature in A is compared to B
          in search of overlaps. Use 'stdin' if passing A with a UNIX pipe.
       -g <path>  Specify a genome file or identifier that defines the order
           and size of the sequences.
       -d  Report the depth at each genome position with 1-based coordinates.
     --dz  Report the depth at each genome position with 0-based coordinates.
           Unlike, -d, this reports only non-zero positions.
     --bg  Report depth as a GRanges.
    --bga  Report depth as a GRanges, as above (i.e., -bg).
           However with this option, regions with zero coverage are also
           reported.
  --split  Treat \"split\" BAM or BED12 entries as distinct BED intervals when
           computing coverage. For BAM files, this uses the CIGAR 'N' and 'D'
           operations to infer the blocks for computing coverage.
           For BED12 files, this uses the BlockCount, BlockStarts, and
           BlockEnds fields (i.e., columns 10,11,12).
 --strand <strand>  Calculate coverage of intervals from a specific strand.
       -5  Calculate coverage of 5' positions (instead of entire interval).
       -3  Calculate coverage of 3' positions (instead of entire interval).
    --max <max>  Combine all positions with a depth >= max into a single bin in
           the histogram.
  --scale <scale>  Scale the coverage by a constant factor.
           Each coverage value is multiplied by this factor before being
           reported. Useful for normalizing coverage by, e.g.,
           reads per million (RPM). [default: 1.0] i.e., unscaled.
     --pc  Calculates coverage of intervals from left point of a pair reads
           to the right point. Works for BAM files only
     --fs <size>  Forces to use fragment size instead of read length.
           Works for BAM files only"

do_bedtools_genomecov <- make_do(R_bedtools_genomecov)
