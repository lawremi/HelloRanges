\name{bedtools_slop}

\alias{bedtools_slop}
\alias{R_bedtools_slop}
\alias{do_bedtools_slop}

\title{bedtools_slop}

\description{
  Widen ranges on the left and/or right side.
}
\usage{
bedtools_slop(cmd = "--help")
R_bedtools_slop(i, b = 0, l = 0, r = 0, s = FALSE, pct = FALSE,
                g = NULL, header = FALSE)
do_bedtools_slop(i, b = 0, l = 0, r = 0, s = FALSE, pct = FALSE,
                 g = NULL, header = FALSE)
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
  \item{b}{
    Widen the same number base pairs in each direction.
  }
  \item{l}{
    The number of base pairs to subtract from the start coordinate.
  }
  \item{r}{
    The number of base pairs to add to the end coordinate.
  }
  \item{s}{
    Define \code{l} and \code{r} based on strand. For example. if
    \code{TRUE}, \code{l=500} for a negative-stranded feature will
    add 500 bp to the end coordinate.
  }
  \item{pct}{
    Define \code{l} and \code{r} as a fraction of the feature
    length. E.g. if used on a 1000bp feature, and \code{l} is 0.50, will
    add 500 bp upstream.
  }
  \item{g}{
    Genome file, identifier or Seqinfo object that defines the order and
    size of the sequences.
  }
  \item{header}{
    Ignored.
  }
}
  
\details{
  As with all commands, there are three interfaces to the
  \command{slop} command:
  \describe{
    \item{\code{bedtools_slop}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_slop}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_slop}}{Evaluates the result of
      \code{R_bedtools_slop}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }

  This is a fairly straight-forward application of
  \code{\link[GenomicRanges]{resize}} and the \code{+} operator on
  GRanges.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  GRanges, or similar, object. In principle, this is an endomorphism.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/slop.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{intra-range-methods} for \code{resize}.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "slop", package="HelloRanges"))
}
## widen on both ends
bedtools_slop("-i a.bed -b 5 -g tiny.genome")
## widen only on the left end
bedtools_slop("-i a.bed -l 5 -r 0 -g tiny.genome")
}
