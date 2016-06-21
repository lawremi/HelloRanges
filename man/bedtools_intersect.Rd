\name{bedtools_intersect}

\alias{bedtools_intersect}
\alias{R_bedtools_intersect}
\alias{do_bedtools_intersect}

\title{bedtools_intersect}

\description{
  Finds and/or counts the intersections between two ranged datasets.
}
\usage{
bedtools_intersect(cmd = "--help")
R_bedtools_intersect(a, b, ubam = FALSE, bed = FALSE, wa = FALSE, wb = FALSE,
                     loj = FALSE, wo = FALSE, wao = FALSE, u = FALSE,
                     c = FALSE, v = FALSE, f = 1e-09, F = 1e-09,
                     r = FALSE, e = FALSE, s = FALSE, S = FALSE,
                     split = FALSE, g = NA_character_, header = FALSE,
                     names = NULL, filenames = FALSE, sortout = FALSE)
do_bedtools_intersect(a, b, ubam = FALSE, bed = FALSE, wa = FALSE, wb = FALSE,
                      loj = FALSE, wo = FALSE, wao = FALSE, u = FALSE,
                      c = FALSE, v = FALSE, f = 1e-09, F = 1e-09,
                      r = FALSE, e = FALSE, s = FALSE, S = FALSE,
                      split = FALSE, g = NA_character_, header = FALSE,
                      names = NULL, filenames = FALSE, sortout = FALSE)
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
    a ranged data structure, such as a GRanges. Each feature in \code{a}
    is compared to \code{b} in search of overlaps. Use \code{"stdin"}
    for input from another process (presumably while running via
    \command{Rscript}). For streaming from a subprocess, prefix the
    command string with \dQuote{<}, e.g.,
    \code{"<grep foo file.bed"}. Any streamed data is assumed to be in
    BED format.
  }
  \item{b}{
    Like \code{a}, except supports multiple datasets, either as a
    vector/list or a comma-separated string. Also supports file glob
    patterns, i.e., strings containing the wildcard, \dQuote{*}.
  }
  \item{ubam}{
    \strong{Not supported yet.} Write uncompressed BAM output. The
    default is write compressed BAM output.
  }
  \item{bed}{
    When using BAM input, return a GRanges (with a \dQuote{blocks}
    column) instead of a GAlignments. For VCF input, return a VRanges
    instead of a VCF object.
  }
  \item{wa}{
    Return the original entry in \code{a} for each overlap.
  }
  \item{wb}{
    Return the original entry in \code{b} for each overlap.
  }
  \item{loj}{
    Perform a \sQuote{left outer join}. That is, for each feature in
    \code{a} report each overlap with \code{b}. If no overlaps are found,
    report an empty range for \code{b} on the \dQuote{.}
    sequence. Implies \code{wa=TRUE} and \code{wb=TRUE}.
  }
  \item{wo}{
    Return the number of base pairs of
    overlap between the two features as the \dQuote{overlap_width}
    metadata column. Implies \code{wa=TRUE} and \code{wb=TRUE}.
  }
  \item{wao}{
    Like \code{wo}, except it additionally implies \code{loj=TRUE}.
  }
  \item{u}{
    Like \code{wa}, except only the unique entries in \code{a} are returned.
  }
  \item{c}{
    Like \code{wa}, except also count the number of hits in \code{b} for
    each range in \code{a} and return the count as the
    \dQuote{overlap_count} metadata column.
  }
  \item{v}{
    Like \code{wa}, except only report those entries in \code{a} that
    have \emph{no} overlap in \code{b}.
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
    Require same strandedness. That is, find the intersect feature in
    \code{b} that overlaps \code{a} on the \emph{same} strand. By default,
    overlaps are reported without respect to strand. Note that this is
    the exact opposite of Bioconductor behavior.
  }
  \item{S}{
    Require opposite strandedness. That is, find the intersect feature in
    \code{b} that overlaps \code{a} on the \emph{opposite} strand. By
    default, overlaps are reported without respect to strand.
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
  \item{header}{
    Ignored.
  }
  \item{names}{
    When using multiple databases, provide an alias for each to use
    instead of their integer index. If a single string, can be
    comma-separated.
  }
  \item{filenames}{
    When using multiple databases, use their complete filename instead
    of their integer index.
  }
  \item{sortout}{
    Sort the result by genomic coordinate.
  }
}
\details{
  As with all commands, there are three interfaces to the
  \command{intersect} command:
  \describe{
    \item{\code{bedtools_intersect}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_intersect}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_intersect}}{Evaluates the result of
      \code{R_bedtools_intersect}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }
  
  This is by far the most complex \command{bedtools} command, and it
  offers a dizzying list of parameters, many of which are redundant or
  mutually exclusive. The complexity of the generated code is highest
  when using the fractional restriction feature, since no such support
  exists in the GenomicRanges overlap routines.  
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  ranges object, the exact type of which depends on the arguments. If
  both \code{wa} and \code{wb} are \code{TRUE}, return a Pairs object
  with the original, matched up ranges, possibly with metadata
  columns. By default, the return value is a GAlignments for BAM input,
  a VCF object for VCF input, or a GRanges for any other type of input.
  If \code{bed} is \code{TRUE}, BAM input is converted to a GRanges,
  containing a \dQuote{blocks} column (encoding the junctions) if the
  input is BAM. If the input is VCF, \code{bed=TRUE} converts the input
  to a VRanges.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{setops-methods} for set operations including
  intersect, \link[GenomicRanges]{findOverlaps-methods} for different
  ways to detect overlaps.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "intersect", package="HelloRanges"))
}

## return intersecting ranges
bedtools_intersect("-a a.bed -b a.bed")
## add count
bedtools_intersect("-a a.bed -b b.bed -c")
## restrict by strand and fraction of overlap
bedtools_intersect("-a a.bed -b b.bed -c -s -f 0.1")
## return original 'a' ranges
bedtools_intersect("-a a.bed -b b.bed -wa")
## return both 'a' and 'b' ranges, along with overlap widths
bedtools_intersect("-a a.bed -b b.bed -wo")
## same as above, except left outer join
bedtools_intersect("-a a.bed -b b.bed -wao")
## consider read junction structure
bedtools_intersect("-a three_blocks.bam -b three_blocks_nomatch.bed -split")
}
