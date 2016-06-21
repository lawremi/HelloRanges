\name{bedtools_jaccard}

\alias{bedtools_jaccard}
\alias{R_bedtools_jaccard}
\alias{do_bedtools_jaccard}

\title{bedtools_jaccard}

\description{
  Compare two sets of genomic regions using the Jaccard statistic,
  defined as the total width of the intersection, divided by the total
  width of the union.
}
\usage{
bedtools_jaccard(cmd = "--help")
R_bedtools_jaccard(a, b, f = 1e-09, F = 1e-09, r = FALSE, e = FALSE,
                   s = FALSE, S = FALSE, split = FALSE)
do_bedtools_jaccard(a, b, f = 1e-09, F = 1e-09, r = FALSE, e = FALSE,
                    s = FALSE, S = FALSE, split = FALSE)
}

\arguments{
  \item{cmd}{
    String of bedtools command line arguments, as they would be entered
    at the shell. There are a few incompatibilities between the
    \CRANpkg{docopt} parser and the bedtools style. See
    \link[=argparsing]{argument parsing}.
  }
  \item{a}{
    Path to a BAM/BED/GFF/VCF/etc file, a BED stream, a file object, or
    a ranged data structure, such as a GRanges. Use \code{"stdin"} for
    input from another process (presumably while running via
    \command{Rscript}). For streaming from a subprocess, prefix the
    command string with \dQuote{<}, e.g.,
    \code{"<grep foo file.bed"}. Any streamed data is assumed to be in
    BED format.
  }
  \item{b}{
    Like \code{a}, except supports multiple datasets, either as a
    vector/list or a comma-separated string. Also supports file glob
    patterns, i.e., strings containing the wildcard, \dQuote{*}.
  }
  \item{f}{
    Minimum overlap required as a fraction of \code{a} [default: any overlap].
  }
  \item{F}{
    Minimum overlap required as a fraction of \code{b} [default: any overlap].
  }
  \item{r}{
    Require that the fraction of overlap be reciprocal for \code{a} and
    \code{b}. In other words, if \code{f} is 0.90 and \code{r} is
    \code{TRUE}, this requires that \code{b} overlap at least 90\% of
    \code{a} and that \code{a} also overlaps at least 90\% of \code{b}.
  }
  \item{e}{
    Require that the minimum fraction be satisfied for \code{a}
    \emph{OR} \code{b}. In other words, if \code{e} is \code{TRUE} with
    \code{f=0.90} and \code{F=0.10} this requires that either 90\% of
    \code{a} is covered OR 10\% of \code{b} is covered. If
    \code{FALSE}, both fractions would have to be satisfied.
  }
  \item{s}{
    Require same strandedness. That is, find the jaccard feature in
    \code{b} that overlaps \code{a} on the \emph{same} strand. By default,
    overlaps are reported without respect to strand. Note that this is
    the exact opposite of Bioconductor behavior.
  }
  \item{S}{
    Require opposite strandedness. That is, find the jaccard feature in
    \code{b} that overlaps \code{a} on the \emph{opposite} strand. By
    default, overlaps are reported without respect to strand.
  }
  \item{split}{
    Treat split BAM (i.e., having an \sQuote{N} CIGAR operation) or
    BED12 entries as compound ranges with gaps, i.e., as GRangesList
    objects.
  }
}
\details{
  As with all commands, there are three interfaces to the
  \command{jaccard} command:
  \describe{
    \item{\code{bedtools_jaccard}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_jaccard}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_jaccard}}{Evaluates the result of
      \code{R_bedtools_jaccard}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }
  
  This is mostly just \code{\link[GenomicRanges]{intersect}} and
  \code{\link[GenomicRanges]{union}}, except when fractional overlap
  restrictions are involved.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  a DataFrame with four columns:
  \item{intersection}{total width of intersection}
  \item{union}{total width of union}
  \item{jaccard}{the jaccard statistic}
  \item{n_intersections}{the number of ranges representing the intersection}
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{setops-methods} for set operations including
  intersect and union.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "jaccard", package="HelloRanges"))
}

## basic
bedtools_jaccard("-a a.bed -b a.bed")
## excluding the gaps in compound ranges
bedtools_jaccard("-a three_blocks_match.bed -b e.bed -split")
## strand and fractional overlap restriction
bedtools_jaccard("-a aMixedStrands.bed -b bMixedStrands.bed -s -f 0.8")
}
