\name{bedtools_coverage}

\alias{bedtools_coverage}
\alias{R_bedtools_coverage}
\alias{do_bedtools_coverage}

\title{bedtools_coverage}

\description{
  Compute the coverage of one or more datasets over a set of query ranges.
}
\usage{
bedtools_coverage(cmd = "--help")
R_bedtools_coverage(a, b, hist = FALSE, d = FALSE, counts = FALSE,
                    f = 1e-09, F = 1e-09, r = FALSE, e = FALSE,
                    s = FALSE, S = FALSE, split = FALSE, g = NA_character_,
                    header = FALSE, sortout = FALSE)
do_bedtools_coverage(a, b, hist = FALSE, d = FALSE, counts = FALSE,
                     f = 1e-09, F = 1e-09, r = FALSE, e = FALSE,
                     s = FALSE, S = FALSE, split = FALSE, g = NA_character_,
                     header = FALSE, sortout = FALSE)
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
    a ranged data structure, such as a GRanges. The coverage is computed
    over these ranges. Use \code{"stdin"} for input from another process
    (presumably while running via \command{Rscript}). For streaming from
    a subprocess, prefix the command string with \dQuote{<}, e.g.,
    \code{"<grep foo file.bed"}. Any streamed data is assumed to be in
    BED format.
  }
  \item{b}{
    Like \code{a}, except supports multiple datasets, either as a
    vector/list or a comma-separated string. Also supports file glob
    patterns, i.e., strings containing the wildcard, \dQuote{*}. The
    coverage is computed by counting how many of these ranges overlap
    positions in \code{a}.
  }
  \item{hist}{
    Report a histogram of coverage for each feature in \code{a} as well as a
    summary histogram for \emph{all} features in \code{a}. See below for
    the structure of the returned table.
  }
  \item{d}{
    Report the depth at each position in each \code{a} feature.
    Positions reported are one based. Each position and depth
    follow the complete \code{a} feature.
  }
  \item{counts}{
    Only report the count of overlaps, not fraction, etc.
    Restricted by \code{f} and \code{r}.
  }
  \item{f}{
    Minimum overlap required as a fraction of A [default: any overlap].
  }
  \item{F}{
    Minimum overlap required as a fraction of B [default: any overlap].
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
    Force strandedness. That is, only count ranges in \code{b} that
    overlap \code{a} on the same strand. By default, coverage is
    computed without respect to strand. Note that this is the exact
    opposite of Bioconductor behavior.
  }
  \item{S}{
    Require opposite strandedness. That is, count the features in
    \code{b} that overlap \code{a} on the \emph{opposite} strand. By
    default, coverage is computed without respect to strand.
  }
  \item{split}{
    Treat split BAM (i.e., having an \sQuote{N} CIGAR operation) or
    BED12 entries as compound ranges with gaps, i.e., as GRangesList
    objects.
  }
  \item{g}{
    A genome file, identifier or Seqinfo object that defines the order
    and size of the sequences.
  }
  \item{header}{
    Ignored.
  }
  \item{sortout}{
    Sort the result by position.
  }
}
  
\details{
  As with all commands, there are three interfaces to the
  \command{coverage} command:
  \describe{
    \item{\code{bedtools_coverage}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_coverage}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_coverage}}{Evaluates the result of
      \code{R_bedtools_coverage}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }

  Typically, we compute coverage with
  \code{\link[GenomicRanges]{coverage}}, but features like fractional
  overlap restriction and histograms add (educational) complexity. One
  key trick is the \code{[,List,GenomicRanges} method, which lets us
  extract coverage vectors for specific regions (see the generated
  code).
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  GRanges object with coverage information. The exact type of
  information depends on the mode:
  \item{default}{Three metadata columns: \dQuote{count} (the number of
    overlapping ranges), \dQuote{covered} (the number of bases covered
    in the query), \dQuote{fraction} (the fraction of bases covered).
  }
  \item{d}{Metadata column \dQuote{coverage} is an RleList with
    position-level coverage (depth). This is what we typically refer to
    as coverage in Bioconductor.
  }
  \item{hist}{
    Metadata column \dQuote{coverage} is a list of DataFrames. Each
    DataFrame contains a histogram of the coverage depth over that
    range, with columns \dQuote{coverage} (the coverage value),
    \dQuote{count} (the number of positions with that coverage),
    \dQuote{len} (the length of the region, all the same) and
    \dQuote{fraction} (the fraction of positions at that coverage).
    There is also a \dQuote{coverage} component on \code{metadata(ans)}
    with the same histogram aggregated over all query ranges.
  }
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{coverage-methods} for ways to compute coverage.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "coverage", package="HelloRanges"))
}

## default behavior
bedtools_coverage("-a a.bed -b b.bed")
## histogram
bedtools_coverage("-a a.bed -b b.bed -hist -g test.genome")
## per-position depth
bedtools_coverage("-a a.bed -b b.bed -d -g test.genome")
}
