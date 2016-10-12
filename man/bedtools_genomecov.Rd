\name{bedtools_genomecov}

\alias{bedtools_genomecov}
\alias{R_bedtools_genomecov}
\alias{do_bedtools_genomecov}

\title{bedtools_genomecov}

\description{
  Compute coverage over the genome. By default, this computes a
  per-chromosome histogram of the coverage, but options allow for
  per-position coverage to be returned in different ways.
}
\usage{
bedtools_genomecov(cmd = "--help")
R_bedtools_genomecov(i, g = NA_character_, d = FALSE, dz = FALSE,
                     bg = FALSE, bga = FALSE, split = FALSE,
                     strand = c("any", "+", "-"), `5` = FALSE, `3` = FALSE,
                     max = NULL, scale = 1, pc = FALSE, fs = NULL)
do_bedtools_genomecov(i, g = NA_character_, d = FALSE, dz = FALSE,
                      bg = FALSE, bga = FALSE, split = FALSE,
                      strand = c("any", "+", "-"), `5` = FALSE, `3` = FALSE,
                      max = NULL, scale = 1, pc = FALSE, fs = NULL)
}

\arguments{
  \item{cmd}{
    String of bedtools command line arguments, as they would be entered
    at the shell. There are a few incompatibilities between the
    \CRANpkg{docopt} parser and the bedtools style. See
    \link[=argparsing]{argument parsing}.
  }
  \item{i}{
    Path to a BAM/BED/GFF/VCF/etc file, a BED stream, a file object, or
    a ranged data structure, such as a GRanges. Use \code{"stdin"} for
    input from another process (presumably while running via
    \command{Rscript}). For streaming from a subprocess, prefix the
    command string with \dQuote{<}, e.g.,
    \code{"<grep foo file.bed"}. Any streamed data is assumed to be in
    BED format.
  }
  \item{g}{
    Genome file, identifier or Seqinfo object that defines the order and
    size of the sequences.
  }
  \item{d}{
    Report the depth at each genome position. This causes a GPos object
    to be returned.
  }
  \item{dz}{
    Same as \code{d}, except the zero coverage positions are dropped.
  }
  \item{bg}{
    Like \code{d}, except returns a GRanges object, which is useful for
    operating on runs of coverage. Zero coverage runs are dropped.
  }
  \item{bga}{
    Like \code{bg}, except the zero coverage runs are retained.
  }
  \item{split}{
    Treat split BAM (i.e., having an \sQuote{N} CIGAR operation) or
    BED12 entries as compound ranges with gaps, i.e., as GRangesList
    objects.
  }
  \item{strand}{
    Calculate coverage of intervals from a specific strand.
  }
  \item{5}{
    Calculate coverage of 5' positions (instead of entire interval).
  }
  \item{3}{
    Calculate coverage of 3' positions (instead of entire interval).
  }
  \item{max}{
    Combine all positions with a depth >= \code{max} into a single bin in
    the histogram.
  }
  \item{scale}{
    Scale the coverage by a constant factor.
    Each coverage value is multiplied by this factor before being
    reported. Useful for normalizing coverage by, e.g.,
    reads per million (RPM).
  }
  \item{pc}{
    Calculates coverage of intervals from left point of a pair reads
    to the right point. Works for BAM files only.
  }
  \item{fs}{
    Forces to use the given fragment size instead of read length.
  }
}
  
\details{
  As with all commands, there are three interfaces to the
  \command{genomecov} command:
  \describe{
    \item{\code{bedtools_genomecov}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_genomecov}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_genomecov}}{Evaluates the result of
      \code{R_bedtools_genomecov}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }

  We typically compute the coverage with
  \code{\link[GenomicRanges]{coverage}}. Computing the histogram
  requires more work.
}
\value{
  A language object containing the compiled \R{} code, evaluating to an
  object that depends on the mode:
  
  \item{default}{A DataFrame that is a per-chromosome histogram of the
    coverage that includes the whole genome margin. Includes columns
    \dQuote{seqnames} (the chromosome name, or \dQuote{genome}),
    \dQuote{coverage} (the coverage value), \dQuote{count} (the count of
    positions covered at that value), \dQuote{len} (the length of the
    chromosome/genome), \dQuote{fraction} (the fraction of bases covered
    at the value).
  }
  \item{d, dz}{A GPos object with per-position coverage values.}
  \item{bg, bga}{A GRanges object with coverage runs.}
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{intra-range-methods} for \code{genomecov}.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "genomecov", package="HelloRanges"))
}
## get coverage runs as a GRanges
bedtools_genomecov("-i y.bed -bg -g test.genome")
## get coverage depth as a GPos, dropping zero values, ignore junctions
bedtools_genomecov("-i three_blocks.bam -dz -split")
## custom fragment size
bedtools_genomecov("-i chip.bam -bg -fs 100")
}
