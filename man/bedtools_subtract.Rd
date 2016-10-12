\name{bedtools_subtract}

\alias{bedtools_subtract}
\alias{R_bedtools_subtract}
\alias{do_bedtools_subtract}

\title{bedtools_subtract}

\description{
  Subtracts one set of ranges from another, either by position or range.
}
\usage{
bedtools_subtract(cmd = "--help")
R_bedtools_subtract(a, b, f = 1e-09, F = 1e-09, r = FALSE, e = FALSE,
                    s = FALSE, S = FALSE, A = FALSE, N = FALSE,
                    g = NA_character_)
do_bedtools_subtract(a, b, f = 1e-09, F = 1e-09, r = FALSE, e = FALSE,
                     s = FALSE, S = FALSE, A = FALSE, N = FALSE,
                     g = NA_character_)
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
    a ranged data structure, such as a GRanges. Each feature in \code{a}
    is compared to \code{b} in search of overlaps. Use \code{"stdin"}
    for input from another process (presumably while running via
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
    Require same strandedness. That is, find the subtract feature in
    \code{b} that overlaps \code{a} on the \emph{same} strand. By default,
    overlaps are reported without respect to strand. Note that this is
    the exact opposite of Bioconductor behavior.
  }
  \item{S}{
    Require opposite strandedness. That is, find the subtract feature in
    \code{b} that overlaps \code{a} on the \emph{opposite} strand. By
    default, overlaps are reported without respect to strand.
  }
  \item{A}{
    Remove entire feature if any overlap. If a feature in \code{a}
    overlaps one in \code{b}, the entire feature is removed.
  }
  \item{N}{
    Same as \code{A=TRUE} except when considering \code{f} the numerator
    in the fraction is the sum of the overlap for all overlapping
    features in \code{b}.
  }
  \item{g}{
    A genome file, identifier or Seqinfo object that defines the order
    and size of the sequences.
  }
}
\details{
  As with all commands, there are three interfaces to the
  \command{subtract} command:
  \describe{
    \item{\code{bedtools_subtract}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_subtract}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_subtract}}{Evaluates the result of
      \code{R_bedtools_subtract}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }

  We typically subtract sets of ranges using
  \code{\link[GenomicRanges]{setdiff}}; however, that will not work
  here, because we cannot merge the ranges in \code{a}.

  The algorithm has two modes: by position (where ranges are clipped)
  and by range (where ranges are discarded entirely). The position mode
  is the default. We find overlaps, optionally restrict them, and for
  each range in \code{a}, we subtract all of the qualifying
  intersections in \code{b}.

  When \code{A} or \code{N} are \code{TRUE}, we use the second mode. In
  the simplest case, that is just
  \code{\link[GenomicRanges]{subsetByOverlaps}} with \code{invert=TRUE},
  but fractional overlap restrictions and \code{N} make that more
  complicated.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  GRanges object, except when \code{A} or \code{N} are \code{TRUE},
  where the value might be a GRanges, GAlignments or VCF object,
  depending on the input.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{setops-methods} for set operations including
  setdiff, \link[GenomicRanges]{findOverlaps-methods} for different
  ways to detect overlaps.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "subtract", package="HelloRanges"))
}

## simple case, position-wise subtraction
bedtools_subtract("-a a.bed -b b.bed")
## fractional overlap restriction
bedtools_subtract("-a a.bed -b b.bed -f 0.5")
## range-wise subtraction
bedtools_subtract("-a a.bed -b b.bed -A -f 0.5")
}
