\name{bedtools_groupby}

\alias{bedtools_groupby}
\alias{R_bedtools_groupby}
\alias{do_bedtools_groupby}

\title{bedtools_groupby}

\description{
    Query sequence from a FASTA file given a set of ranges, including
    compound regions like transcripts and junction reads. This assumes
    the sequence is DNA.
}
\usage{
    bedtools_groupby(cmd = "--help")
    R_bedtools_groupby(i, g = 1:3, c, o = "sum", delim=",")
    do_bedtools_groupby(i, g = 1:3, c, o = "sum", delim=",")
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
    Column index(es) for grouping the input. Columns may be
    comma-separated. By default, the grouping is by range.
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
    \code{o}. See the details for the available operations.
  }
  \item{delim}{
    Delimiter character used to collapse strings.
  }
}

\details{
  As with all commands, there are three interfaces to the
  \command{groupby} command:
  \describe{
    \item{\code{bedtools_groupby}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_groupby}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_groupby}}{Evaluates the result of
      \code{R_bedtools_groupby}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }
  
  The workhorse for aggregation in \R{} is
  \code{\link[stats]{aggregate}} and we have extended its
  interface to make it more convenient. See
  \code{\link[S4Vectors]{aggregate}} for details.

  The following operations are supported (with \R{} translation):
  \describe{
    \item{sum}{\code{sum(X)}}
    \item{min}{\code{min(X)}}
    \item{max}{\code{max(X)}}
    \item{absmin}{\code{min(abs(X))}}
    \item{absmax}{\code{max(abs(X))}}
    \item{mean}{\code{mean(X)}}
    \item{median}{\code{median(X)}}
    \item{mode}{\code{\link{distmode}(X)}}
    \item{antimode}{\code{distmode(X, anti=TRUE)}}
    \item{collapse}{\code{unstrsplit(X, delim)}}
    \item{distinct}{\code{unstrsplit(unique(X), delim)}}
    \item{count}{\code{lengths(X)}}
    \item{count_distinct}{\code{lengths(unique(X))}}
    \item{sstdev}{\code{sd(X)}}
    \code{freq}{\code{table(X)}}
    \code{first}{\code{drop(phead(X, 1L))}}
    \code{last}{\code{drop(ptail(X, 1L))}}
  }

  For the sake of simplicity, and because the use cases are not clear,
  we do not support aggregation of every column. Here are some of the
  restrictions:
  \itemize{
    \item{
      No support for the last column of GFF (the ragged list of
      attributes).
    }
    \item{
      No support for the INFO, FORMAT and GENO fields of VCF.
    }
    \item{
      No support for the FLAG field of BAM (\command{bedtools} does
      not support this either).
    }
  }
  
}
\note{
  We admit that using column subscripts for \code{c} makes code hard
  to read. All the more reason to just write \R{} code.  
}
\value{
    A language object containing the compiled \R{} code, generally
    evaluating to a DataFrame, with a column for each grouping variable
    and each summarized variable. As a special case, if there are no
    grouping variables specified, then the grouping is by range, and an
    aggregated GRanges is returned.
}
\references{
    \url{http://bedtools.readthedocs.io/en/latest/content/tools/groupby.html}
}
\author{
    Michael Lawrence
}

\seealso{
    \code{\link[S4Vectors]{aggregate-methods}} for general aggregation.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "groupby", package="HelloRanges"))
}
    ## aggregation by range
    bedtools_groupby("-i values3.header.bed -c 5")
    ## average variant qualities by chromosome and reference base
\dontrun{
    indexTabix(bgzip("a_vcfSVtest.vcf", overwrite=TRUE), "vcf")
}
    bedtools_groupby("-i a_vcfSVtest.vcf.bgz -g 1,4 -c 6 -o mean")
}
