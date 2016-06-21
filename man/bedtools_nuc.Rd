\name{bedtools_nuc}

\alias{bedtools_nuc}
\alias{R_bedtools_nuc}
\alias{do_bedtools_nuc}

\title{bedtools_nuc}

\description{
    Summarize DNA sequences over the specified ranges.
}
\usage{
    bedtools_nuc(cmd = "--help")
    R_bedtools_nuc(fi, bed, s = FALSE, pattern = NULL, fullHeader = FALSE)
    do_bedtools_nuc(fi, bed, s = FALSE, pattern = NULL, fullHeader = FALSE)
}

\arguments{
    \item{cmd}{
      String of bedtools command line arguments, as they would be entered
      at the shell. There are a few incompatibilities between the
      \CRANpkg{docopt} parser and the bedtools style. See
      \link[=argparsing]{argument parsing}.
    }
    \item{fi}{
        Path to a FASTA file, or an XStringSet.
    }
    \item{bed}{
      Path to a BAM/BED/GFF/VCF/etc file, a BED stream, a file object, or
      a ranged data structure, such as a GRanges, as the query. Use
      \code{"stdin"} for input from another process (presumably while
      running via \command{Rscript}). For streaming from a subprocess,
      prefix the command string with \dQuote{<}, e.g.,
      \code{"<grep foo file.bed"}. Any streamed data is assumed to be in
      BED format.
    }
    \item{s}{
      Force strandedness. If the feature occupies the antisense strand,
      the sequence will be reverse complemented.
    }
    \item{pattern}{
      Optional sequence pattern to count in each subsequence.
    }
    \item{fullHeader}{
      Use the full FASTA header as the names. By default, use just the
      first word.
    }
}
\details{
    As with all commands, there are three interfaces to the
    \command{nuc} command:
    \describe{
        \item{\code{bedtools_nuc}}{Parses the bedtools command line and
            compiles it to the equivalent \R{} code.
        }
        \item{\code{R_bedtools_nuc}}{Accepts \R{} arguments
            corresponding to the command line arguments and compiles the
            equivalent \R{} code.
        }
        \item{\code{do_bedtools_nuc}}{Evaluates the result of
            \code{R_bedtools_nuc}. Recommended \strong{only} for
            demonstration and testing. It is best to integrate the compiled
            code into an R script, after studying it.
        }
    }
    
    Computes AT/GC percentage and counts each type of base. Relies on
    Biostrings utilities like \code{\link[Biostrings]{letterFrequency}}
    and \code{\link[Biostrings]{alphabetFrequency}}. The counting of
    \code{pattern} occurrences uses
    \code{\link[Biostrings]{vcountPattern}}.
}
\value{
    A language object containing the compiled \R{} code, evaluating to a
    DataFrame with summary statistics including the AC and GT
    percentage, and the counts of each type of base. Also includes the
    count of \code{pattern}, if specified.
}
\references{
    \url{http://bedtools.readthedocs.io/en/latest/content/tools/nuc.html}
}
\author{
    Michael Lawrence
}

\seealso{
  \link[Biostrings]{letterFrequency} for summarizing sequences,
  \link[Biostrings]{matchPattern} for pattern matching.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "nuc", package="HelloRanges"))
}
    ## default behavior, note the two dashes in '--fi'
    bedtools_nuc("--fi test.fasta -bed a.bed")
    ## with pattern counting
    bedtools_nuc("--fi test.fasta -bed a.bed -pattern ATA")
}
