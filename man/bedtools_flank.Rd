\name{bedtools_flank}

\alias{bedtools_flank}
\alias{R_bedtools_flank}
\alias{do_bedtools_flank}

\title{bedtools_flank}

\description{
  Compute flanking regions.
}
\usage{
bedtools_flank(cmd = "--help")
R_bedtools_flank(i, b = 0, l = 0, r = 0, s = FALSE, pct = FALSE,
                 g = NULL, header = FALSE)
do_bedtools_flank(i, b = 0, l = 0, r = 0, s = FALSE, pct = FALSE,
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
    Increase the BED/GFF/VCF range by the same number base pairs in each
    direction. Integer.
  }
  \item{l}{
    The number of base pairs to subtract from the start coordinate. Integer.
  }
  \item{r}{
    The number of base pairs to add to the end coordinate. Integer.
  }
  \item{s}{
    Define \code{l} and \code{r} based on strand. For example. if used,
    \code{l} is 500 for a negative-stranded feature, it will add 500 bp
    to the end coordinate.
  }
  \item{pct}{
    Define \code{l} and \code{r} as a fraction of the feature
    length. E.g. if used on a 1000bp feature, and \code{l} is 0.50, will
    add 500 bp upstream..
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
  \command{flank} command:
  \describe{
    \item{\code{bedtools_flank}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_flank}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_flank}}{Evaluates the result of
      \code{R_bedtools_flank}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }

  We compute the flanks with \code{\link[GenomicRanges]{flank}},
  although \code{flank} only computes one side at a time, so we may call
  it multiple times.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  GRanges object.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/flank.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{intra-range-methods} for \code{flank}.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "flank", package="HelloRanges"))
}
## 5 on both sides
r <- bedtools_flank("-i a.bed -b 5 -g tiny.genome")
## 5 on left side
bedtools_flank("-i a.bed -l 5 -r 0 -g tiny.genome")
## define left/right in terms of transcription direction
bedtools_flank("-i a.bed -l 5 -r 0 -s -g tiny.genome")
}
