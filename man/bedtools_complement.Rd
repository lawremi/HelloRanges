\name{bedtools_complement}

\alias{bedtools_complement}
\alias{R_bedtools_complement}
\alias{do_bedtools_complement}

\title{bedtools_complement}

\description{
  Finds regions of the genome that are not covered by a genomic dataset.
}
\usage{
bedtools_complement(cmd = "--help")
R_bedtools_complement(i, g)
do_bedtools_complement(i, g)
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
     A genome file, identifier or Seqinfo object that defines the order
     and size of the sequences.
  }
}
  
\details{
  As with all commands, there are three interfaces to the
  \command{complement} command:
  \describe{
    \item{\code{bedtools_complement}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_complement}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_complement}}{Evaluates the result of
      \code{R_bedtools_complement}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }

  The generated code is subtracts, via
  \code{\link[GenomicRanges]{setdiff}}, the ranges from the set of
  ranges representing the entire genome.

  While it may be tempting to call \code{\link[GenomicRanges]{gaps}}
  instead, it is very unlikely to behave as expected. The GenomicRanges
  set operations treat all three strand values (+, -, *) as separate
  spaces. \code{gaps} takes as its universe the genome on all three
  strands, rather than just the \dQuote{*} strand, resulting in
  extraneous stranded ranges.  
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  GRanges object with the complementary ranges.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/complement.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{setops-methods} for the various set operations.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "coverage", package="HelloRanges"))
}
bedtools_complement("-i a.bed -g test.genome")
}
