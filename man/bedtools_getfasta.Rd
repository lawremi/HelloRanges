\name{bedtools_getfasta}

\alias{bedtools_getfasta}
\alias{R_bedtools_getfasta}
\alias{do_bedtools_getfasta}

\title{bedtools_getfasta}

\description{
    Query sequence from a FASTA file given a set of ranges, including
    compound regions like transcripts and junction reads. This assumes
    the sequence is DNA.
}
\usage{
    bedtools_getfasta(cmd = "--help")
    R_bedtools_getfasta(fi, bed, s = FALSE, split = FALSE)
    do_bedtools_getfasta(fi, bed, s = FALSE, split = FALSE)
}

\arguments{
  \item{cmd}{
    String of bedtools command line arguments, as they would be entered
    at the shell. There are a few incompatibilities between the
    \CRANpkg{docopt} parser and the bedtools style. See
    \link[=argparsing]{argument parsing}.
  }
  \item{fi}{
    Path to a FASTA file, or an XStringSet object.
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
  \item{split}{
    Given BED12 or BAM input, extract and concatenate the sequences from
    the blocks (e.g., exons).
  }
}

\details{
    As with all commands, there are three interfaces to the
    \command{getfasta} command:
    \describe{
        \item{\code{bedtools_getfasta}}{Parses the bedtools command line and
            compiles it to the equivalent \R{} code.
        }
        \item{\code{R_bedtools_getfasta}}{Accepts \R{} arguments
            corresponding to the command line arguments and compiles the
            equivalent \R{} code.
        }
        \item{\code{do_bedtools_getfasta}}{Evaluates the result of
            \code{R_bedtools_getfasta}. Recommended \strong{only} for
            demonstration and testing. It is best to integrate the compiled
            code into an R script, after studying it.
        }
    }

    It is recommended to retrieve reference sequence using a
    \pkg{BSgenome} package, either custom or provided by
    Bioconductor. Call \code{\link[Biostrings]{getSeq}} to query for
    specific regions of the BSgenome object. If one must access a file,
    consider converting it to 2bit or FA (razip) format for indexed
    access using \code{\link[BiocIO]{import}} and its \code{which}
    argument.

    But if one must access a FASTA file, we need to read all of it with
    \code{\link[Biostrings]{readDNAStringSet}} and extract regions using
    \code{x[gr]}, where \code{gr} is a GRanges or GRangesList.
}
\value{
    A language object containing the compiled \R{} code, evaluating to a
    DNAStringSet object.
}
\references{
    \url{http://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html}
}
\author{
    Michael Lawrence
}

\seealso{
    \link[Biostrings]{getSeq}, the primary sequence query interface.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "getfasta", package="HelloRanges"))
}
    ## simple query
    bedtools_getfasta("--fi t.fa -bed blocks.bed")
    ## get spliced transcript/read sequence
    bedtools_getfasta("--fi t.fa -bed blocks.bed -split")
}
