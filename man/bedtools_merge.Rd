\name{bedtools_merge}

\alias{bedtools_merge}
\alias{R_bedtools_merge}
\alias{do_bedtools_merge}

\title{bedtools_merge}

\description{
  Collapse overlapping and adjacent ranges into a single range, i.e.,
  \code{\link{reduce}} the ranges. Then, group the original ranges by
  reduced range and aggregate. By default, the scores are summed.
}
\usage{
bedtools_merge(cmd = "--help")
R_bedtools_merge(i, s = FALSE, S = c("any", "+", "-"), d = 0L, c = NULL,
                 o = "sum", delim = ",")
do_bedtools_merge(i, s = FALSE, S = c("any", "+", "-"), d = 0L, c = NULL,
                  o = "sum", delim = ",")
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
    BED format. These are the ranges that are merged.
  }
  \item{s}{
    Require same strandedness. That is, find the jaccard feature in
    \code{b} that overlaps \code{a} on the \emph{same} strand. By default,
    overlaps are reported without respect to strand. Note that this is
    the exact opposite of Bioconductor behavior.
  }
  \item{S}{
    Force merge for one specific strand only.  Follow with + or - to
    force merge from only the forward or reverse strand, respectively.
    By default, merging is done without respect to strand.
  }
  \item{d}{
    Maximum distance between features allowed for features to be
    merged. Default is 0. That is, overlapping and/or book-ended
    features are merged.
  }
  \item{c}{
    Specify columns (by integer index) from the input file to operate
    upon (see \code{o} option, below).  Multiple columns can be
    specified in a comma-delimited list.
  }
  \item{o}{
    Specify the operations (by name) that should be applied to the
    columns indicated in \code{c}. Multiple operations can be specified
    in a comma-delimited list. Recycling is used to align \code{c} and
    \code{o}. See \code{\link{bedtools_groupby}} for the available
    operations. Defaults to the \dQuote{sum} operation.
  }
  \item{delim}{
    Delimiter character used to collapse strings.
  }
}
\details{
  As with all commands, there are three interfaces to the
  \command{merge} command:
  \describe{
    \item{\code{bedtools_merge}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_merge}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_merge}}{Evaluates the result of
      \code{R_bedtools_merge}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }
  
  The workhorse for reduction is
  \code{\link[GenomicRanges]{reduce}}. Passing \code{with.revmap=TRUE}
  to \code{reduce} causes it to return a list of integers, which can be
  passed directly to \code{\link[S4Vectors]{aggregate}} to aggregate the
  original ranges.

  Since the grouping information is preserved in the result, this
  function serves as a proxy for \command{bedtools cluster}.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  DataFrame with a \dQuote{grouping} column corresponding to
  \code{as(hits, "List")}, and a column for each summary.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/merge.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link{bedtools_groupby} for more details on bedtools-style aggregation,
  \link[GenomicRanges]{reduce} for merging,
  \link[S4Vectors]{aggregate-methods} for aggregating.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "merge", package="HelloRanges"))
}
## default behavior, sum the score
bedtools_merge("-i a.bed")
## count the seqnames
bedtools_merge("-i a.bed -c 1 -o count")
## collapse the names using "|" as the delimiter
bedtools_merge("-i a.names.bed -delim \"|\" -c 4  -o collapse")
## collapse the names and sum the scores
bedtools_merge("-i a.full.bed -c 4,5  -o collapse,sum")
## count and sum the scores
bedtools_merge("-i a.full.bed -c 5  -o count,sum")
## only merge the positive strand features
bedtools_merge("-i a.full.bed -S +")
}
