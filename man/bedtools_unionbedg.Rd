\name{bedtools_unionbedg}

\alias{bedtools_unionbedg}
\alias{R_bedtools_unionbedg}
\alias{do_bedtools_unionbedg}

\title{bedtools_unionbedg}

\description{
    Summarize the ranges according to \code{\link[GenomicRanges]{disjoin}}
    and construct a matrix of scores (disjoint range by
    sample/file). Empty cells are filled with \code{NA}.
}
\usage{
    bedtools_unionbedg(cmd = "--help")
    R_bedtools_unionbedg(i, header=FALSE, names=NULL, g=NA_character_,
                          empty=FALSE)
    do_bedtools_unionbedg(i, header=FALSE, names=NULL, g=NA_character_,
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
    \command{unionbedg} command:
    \describe{
        \item{\code{bedtools_unionbedg}}{Parses the bedtools command line and
            compiles it to the equivalent \R{} code.
        }
        \item{\code{R_bedtools_unionbedg}}{Accepts \R{} arguments
            corresponding to the command line arguments and compiles the
            equivalent \R{} code.
        }
        \item{\code{do_bedtools_unionbedg}}{Evaluates the result of
            \code{R_bedtools_unionbedg}. Recommended \strong{only} for
            demonstration and testing. It is best to integrate the compiled
            code into an R script, after studying it.
        }
    }
    
    This is essentially the same operation as
    \code{\link{bedtools_multiinter}}, except we build a score matrix
    and embed it into a SummarizedExperiment. This is a bit tricky and
    relies on the \code{\link{as.matrix,AtomicList-method}} coercion.
}
\value{
    A language object containing the compiled \R{} code, evaluating to a
    RangedSummarizedExperiment with an assay called \dQuote{score}.
}
\references{
    \url{http://bedtools.readthedocs.io/en/latest/content/tools/unionbedg.html}
}
\author{
    Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{disjoin} for forming disjoint ranges,
  \link[SummarizedExperiment]{RangedSummarizedExperiment-class} for
  SummarizedExperiment objects.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "unionbedg", package="HelloRanges"))
}

## combine three samples
bedtools_unionbedg("-i a.bedGraph,b.bedGraph,c.bedGraph -names A,B,C")
## include empty ranges (filled with NAs)
bedtools_unionbedg("-i a.bedGraph,b.bedGraph,c.bedGraph -names A,B,C -empty -g test.genome")
}
