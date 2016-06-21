\name{bedtools_map}

\alias{bedtools_map}
\alias{R_bedtools_map}
\alias{do_bedtools_map}

\title{bedtools_map}

\description{
  Group ranges by overlap with query ranges and aggregate. By default,
  the scores are summed.
}
\usage{
bedtools_map(cmd = "--help")
R_bedtools_map(a, b, c = "5", o = "sum", f = 1e-09, F = 1e-09,
               r = FALSE, e = FALSE, s = FALSE, S = FALSE, header = FALSE,
               split = FALSE, g = NA_character_, delim=",")
do_bedtools_map(a, b, c = "5", o = "sum", f = 1e-09, F = 1e-09,
                r = FALSE, e = FALSE, s = FALSE, S = FALSE, header = FALSE,
                split = FALSE, g = NA_character_, delim=",")
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
    BED format. Windows are generated with each range. Exclusive with
    \code{g}. A summary of \code{b} is computed for each range.
  }
  \item{b}{
    Like \code{a}, except supports multiple datasets, either as a
    vector/list or a comma-separated string. Also supports file glob
    patterns, i.e., strings containing the wildcard, \dQuote{*}. Ranges
    that map to the same range in \code{a} are aggregated.
  }
  \item{c}{
    Specify columns (by integer index) from the input file to operate
    upon (see \code{o} option, below).  Multiple columns can be
    specified in a comma-delimited list. Defaults to the score column.
  }
  \item{o}{
    Specify the operations (by name) that should be applied to the
    columns indicated in \code{c}. Multiple operations can be specified
    in a comma-delimited list. Recycling is used to align \code{c} and
    \code{o}. See \code{\link{bedtools_groupby}} for the available
    operations. Defaults to the \dQuote{sum} operation.
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
  \item{header}{
    Ignored.
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
  \item{delim}{
    Delimiter character used to collapse strings.
  }
}
\details{
  As with all commands, there are three interfaces to the
  \command{map} command:
  \describe{
    \item{\code{bedtools_map}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_map}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_map}}{Evaluates the result of
      \code{R_bedtools_map}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }
  
  Computing overlaps with \code{\link[GenomicRanges]{findOverlaps}}
  generates a \link[S4Vectors]{Hits} object, which we can pass directly
  to \code{\link[S4Vectors]{aggregate}} to aggregate the subject
  features that overlap the same range in the query.

  There are several commands in the \command{bedtools} suite that might
  be approximately implemented by passing multiple files to \code{b} and
  specifying the aggregate expression \code{table(b)}. That counts how
  many ranges from each database/sample overlap a given query. The
  covered commands are: \command{bedtools annotate -counts},
  \code{bedtools multicov} and \code{bedtools tag}.
}
\note{
  We do not support the \command{bedtools} \code{null} argument, because
  it seems more sensible to just let \R{} decide on the value of statistics
  when a group is empty.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  DataFrame with a \dQuote{grouping} column corresponding to
  \code{as(hits, "List")}, and a column for each summary.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/map.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{findOverlaps-methods} for finding hits,
  \link[S4Vectors]{Hits-class} for manipulating them,
  \link[S4Vectors]{aggregate-methods} for aggregating them.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "map", package="HelloRanges"))
}

## default behavior
bedtools_map("-a ivls.bed -b values.bed")
## take the mode of the scores
bedtools_map("-a ivls.bed -b values.bed -o mode")
## collapse the chromosome names
bedtools_map("-a ivls.bed -b test.gff2 -c 1 -o collapse")
## collapse the names, restricted by fractional overlap
bedtools_map("-a ivls2.bed -b values5.bed -c 4 -o collapse -f 0.7")
}
