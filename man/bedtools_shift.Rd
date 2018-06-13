\name{bedtools_shift}

\alias{bedtools_shift}
\alias{R_bedtools_shift}
\alias{do_bedtools_shift}

\title{bedtools_shift}

\description{
  Compute shifting regions.
}
\usage{
bedtools_shift(cmd = "--help")
R_bedtools_shift(i, s = 0, m = 0, p = 0, pct = FALSE, g = NULL, header = FALSE)
do_bedtools_shift(i, s = 0, m = 0, p = 0, pct = FALSE, g = NULL, header = FALSE)
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
  \item{s}{
    Amount to shift all features.
  }
  \item{m}{
    Amount to shift negative strand features.
  }
  \item{p}{
    Amount to shift positive strand features.
  }
  \item{pct}{
    Define \code{l} and \code{r} as a fraction of the feature
    length. E.g. if used on a 1000bp feature, and \code{l} is 0.50, will
    shift 500 bp upstream..
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
  \command{shift} command:
  \describe{
    \item{\code{bedtools_shift}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_shift}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_shift}}{Evaluates the result of
      \code{R_bedtools_shift}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }

  This is a fairly straight-forward application of
  \code{\link[IRanges]{shift}}.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  GRanges, or similar, object. In principle, this is an endomorphism.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/shift.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{intra-range-methods} for \code{shift}.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "shift", package="HelloRanges"))
}
## shift all ranges by 5
bedtools_shift("-i a.bed -s 5 -g tiny.genome")
## shift only the negative strand features by 5
bedtools_shift("-i a.bed -p 0 -m 5 -g tiny.genome")
}
