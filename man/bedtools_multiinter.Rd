\name{bedtools_multiinter}

\alias{bedtools_multiinter}
\alias{R_bedtools_multiinter}
\alias{do_bedtools_multiinter}

\title{bedtools_multiinter}

\description{
  Summarize the ranges according to \code{\link[GenomicRanges]{disjoin}}
  and annotate each disjoint range with the samples that overlap the
  range.
}
\usage{
bedtools_multiinter(cmd = "--help")
R_bedtools_multiinter(i, header=FALSE, names=NULL, g=NA_character_,
                      empty=FALSE)
do_bedtools_multiinter(i, header=FALSE, names=NULL, g=NA_character_,
                       empty=FALSE)
}

\arguments{
  \item{cmd}{
    String of bedtools command line arguments, as they would be entered
    at the shell. There are a few incompatibilities between the
    \CRANpkg{docopt} parser and the bedtools style. See
    \link[=argparsing]{argument parsing}.
  }
  \item{i}{
    Paths to BAM/BED/GFF/VCF/etc files (vector or comma-separated), or a
    list of objects.
  }
  \item{header}{
    Ignored.
  }
  \item{names}{
    Provide an alias for each to use for each \code{i} instead of their
    integer index. If a single string, can be comma-separated.
  }
  \item{g}{
    A genome file, identifier or Seqinfo object that defines the order
    and size of the sequences.
  }
  \item{empty}{
    Report empty regions (i.e., regions not covered in any of the
    files). This essentially yields a partitioning of the genome (and
    thus requires \code{g} to be specified).
  }
}
\details{
  As with all commands, there are three interfaces to the
  \command{multiinter} command:
  \describe{
    \item{\code{bedtools_multiinter}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_multiinter}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_multiinter}}{Evaluates the result of
      \code{R_bedtools_multiinter}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }
  
  The workhorse is
  \code{\link[GenomicRanges]{disjoin}}. Passing \code{with.revmap=TRUE}
  to \code{disjoin} causes it to return a list of integers, which we use
  to extract the sample identifiers. The \code{empty} case requires a
  bit more code, because we have to combine the disjoint ranges with the
  gaps.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  GRanges with a column \dQuote{i} indicating the sample memberships.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/multiinter.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{disjoin} for forming disjoint ranges.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "multiinter", package="HelloRanges"))
}
## default behavior
bedtools_multiinter("-i a.bed,b.bed,c.bed")
## custom names
bedtools_multiinter("-i a.bed,b.bed,c.bed -names A,B,C")
## include empty regions, i.e., partition the genome
bedtools_multiinter("-i a.bed,b.bed,c.bed -names A,B,C -empty -g test.genome")
}
