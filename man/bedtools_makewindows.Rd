\name{bedtools_makewindows}

\alias{bedtools_makewindows}
\alias{R_bedtools_makewindows}
\alias{do_bedtools_makewindows}

\title{bedtools_makewindows}

\description{
  Generate a partitioning/tiling or set of sliding windows over the
  genome or a set of ranges.
}
\usage{
bedtools_makewindows(cmd = "--help")
R_bedtools_makewindows(b, g = NA_character_, w, s, n)
do_bedtools_makewindows(b, g = NA_character_, w, s, n)
}

\arguments{
  \item{cmd}{
    String of bedtools command line arguments, as they would be entered
    at the shell. There are a few incompatibilities between the
    \CRANpkg{docopt} parser and the bedtools style. See
    \link[=argparsing]{argument parsing}.
  }
  \item{b}{
    Path to a BAM/BED/GFF/VCF/etc file, a BED stream, a file object, or
    a ranged data structure, such as a GRanges. Use \code{"stdin"} for
    input from another process (presumably while running via
    \command{Rscript}). For streaming from a subprocess, prefix the
    command string with \dQuote{<}, e.g.,
    \code{"<grep foo file.bed"}. Any streamed data is assumed to be in
    BED format. Windows are generated with each range. Exclusive with
    \code{g}.
  }
  \item{g}{
    A genome file, identifier or Seqinfo object that defines the order
    and size of the sequences. Specifying this generates windows over
    the genome. Exclusive with \code{b}.
  }
  \item{w}{
    Window size, exclusive with \code{n}.
  }
  \item{s}{
    Step size (generates sliding windows).
  }
  \item{n}{
    Number of windows, exclusive with \code{w}.
  }
}
\details{
  As with all commands, there are three interfaces to the
  \command{makewindows} command:
  \describe{
    \item{\code{bedtools_makewindows}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_makewindows}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_makewindows}}{Evaluates the result of
      \code{R_bedtools_makewindows}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }
  
  We view the generation of a partitioning (or tiling) as a distinct use
  case from the generation of sliding windows. The two use cases
  correspond to the \code{\link[IRanges]{tile}} and
  \code{\link[IRanges]{slidingWindows}} functions, respectively.
}
\value{
  A language object containing the compiled \R{} code, evaluating to a a
  GRangesList containing the windows for each range (or chromosome).
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/makewindows.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{tile-methods} for generating windows.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "makewindows", package="HelloRanges"))
}

## tiles of width 5000
bedtools_makewindows("-b input.bed -w 5000")
## sliding windows, 5kb wide, every 2kb
bedtools_makewindows("-b input.bed -w 5000 -s 2000")
## 3 tiles in each range
bedtools_makewindows("-b input.bed -n 3")
## 3 tiles for each chromosome of the genome
bedtools_makewindows("-g test.genome -n 3")
}
