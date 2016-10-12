\name{bedtools_closest}

\alias{bedtools_closest}
\alias{R_bedtools_closest}
\alias{do_bedtools_closest}

\title{bedtools_closest}

\description{
  Finds the features in one dataset that are closest to those in
  another. Supports restriction by strand, upstream, downstream,
  and overlap. There are several methods for resolving ties. Optionally
  returns the distance.
}
\usage{
bedtools_closest(cmd = "--help")
R_bedtools_closest(a, b, s = FALSE, S = FALSE, d = FALSE,
                   D = c("none", "ref", "a", "b"), io = FALSE, iu = FALSE, 
                   id = FALSE, fu = FALSE, fd = FALSE,
                   t = c("all", "first", "last"), mdb = c("each", "all"), 
                   k = 1L, names = NULL, filenames = FALSE, N = FALSE)
do_bedtools_closest(a, b, s = FALSE, S = FALSE, d = FALSE,
                    D = c("none", "ref", "a", "b"), io = FALSE, iu = FALSE, 
                    id = FALSE, fu = FALSE, fd = FALSE,
                    t = c("all", "first", "last"), mdb = c("each", "all"), 
                    k = 1L, names = NULL, filenames = FALSE, N = FALSE)
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
    is compared to \code{b} in search of nearest neighbors. Use \code{"stdin"}
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
  \item{s}{
    Require same strandedness. That is, find the closest feature in
    \code{b} that overlaps \code{a} on the \emph{same} strand. By default,
    overlaps are reported without respect to strand. Note that this is
    the exact opposite of Bioconductor behavior.
  }
  \item{S}{
    Require opposite strandedness. That is, find the closest feature in
    \code{b} that overlaps \code{a} on the \emph{opposite} strand. By
    default, overlaps are reported without respect to strand.
  }
  \item{d}{
    In addition to the closest feature in \code{b}, report its distance
    to \code{a} as an extra column. The reported distance for
    overlapping features will be 0.
  }
  \item{D}{
    Like \code{d}, report the closest feature in \code{b}, and its
    distance to \code{a} as an extra column. However unlike \code{d},
    \code{D} conveys a notion of upstream that is useful with other
    arguments. See details.
  }
  \item{io}{
    Ignore features in \code{b} that overlap \code{a}. That is, we want
    close, yet not touching features only.
  }
  \item{iu}{
    Ignore features in \code{b} that are upstream of features in
    \code{a}. This option requires \code{D} and follows its orientation
    rules for determining what is \dQuote{upstream}.
  }
  \item{id}{
    Ignore features in \code{b} that are downstream of features in
    \code{a}. This option requires \code{D} and follows its orientation
    rules for determining what is \dQuote{downstream}.
  }
  \item{fu}{
    Choose first from features in \code{b} that are upstream of features
    in \code{a}. This option requires \code{D} and follows its
    orientation rules for determining what is \dQuote{upstream}.
  }
  \item{fd}{
    Choose first from features in \code{b} that are downstream of
    features in \code{a}. This option requires \code{D} and follows its
    orientation rules for determining what is \dQuote{downstream}.
  }
  \item{t}{
    Specify how ties for closest feature should be handled. This occurs
    when two features in \code{b} have exactly the same
    \dQuote{closeness} with \code{a}. By default, all such features in
    \code{b} are reported. The modes options are \dQuote{all},
    \dQuote{first} and \dQuote{last}.
  }
  \item{mdb}{
    How multiple databases should be resolved, either \dQuote{each}
    (find closest in each \code{b} dataset independently) or
    \dQuote{all} (combine all \code{b} datasets prior to the search).
  }
  \item{k}{
    \strong{Not supported yet.} Report the \code{k} closest
    hits. Default is 1. If \code{t} is \code{"all"}, all ties will still
    be reported.
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
  \item{N}{
    \strong{Not yet supported}, but related use cases are often solved by
    passing a single argument to
    \code{\link[GenomicRanges]{nearest}}. Require that the query and the
    closest hit have different names. For BED, the 4th column is
    compared.
  }
}
\details{
  As with all commands, there are three interfaces to the
  \command{closest} command:
  \describe{
    \item{\code{bedtools_closest}}{Parses the bedtools command line and
      compiles it to the equivalent \R{} code.
    }
    \item{\code{R_bedtools_closest}}{Accepts \R{} arguments
      corresponding to the command line arguments and compiles the
      equivalent \R{} code.
    }
    \item{\code{do_bedtools_closest}}{Evaluates the result of
      \code{R_bedtools_closest}. Recommended \strong{only} for
      demonstration and testing. It is best to integrate the compiled
      code into an R script, after studying it.
    }
  }
  
  The generated code includes calls to utilities like
  \link[GenomicRanges]{nearest}, \link[GenomicRanges]{precede} and
  \link[GenomicRanges]{follow}. \link[GenomicRanges]{nearest} lacks the
  ability to restrict its search by direction/overlap, so some complex
  code is generated to support all of the argument combinations.

  Arguments \code{io}, \code{iu}, \code{id}, \code{fu}, and \code{fd}
  require a notion of upstream/downstream to be defined via \code{D},
  which accepts one of these values:
  \describe{
    \item{ref}{
      Report distance with respect to the reference genome.  B features
      with a lower (start, stop) are \dQuote{upstream}.
    }
    \item{a}{
      Report distance with respect to A.  When A is on the - strand,
      \dQuote{upstream} means B has a higher (start,stop).
    }
    \item{b}{
      Report distance with respect to B.  When B is on the -
      strand, \dQuote{upstream} means A has a higher (start,stop).
    }
  }
  
}
\value{
  A language object containing the compiled \R{} code, evaluating to a
  Pairs object with the closest hits from \code{a} and \code{b}. If
  \code{d} or \code{D} is \code{TRUE}, has a metadata column called
  \dQuote{distance}.
}
\references{
  \url{http://bedtools.readthedocs.io/en/latest/content/tools/closest.html}
}
\author{
  Michael Lawrence
}

\seealso{
  \link[GenomicRanges]{nearest-methods} for the various ways to find the
  nearest ranges.
}

\examples{
\dontrun{
setwd(system.file("unitTests", "data", "closest", package="HelloRanges"))
}

## basic
bedtools_closest("-a a.bed -b b.bed -d")
## strand-specific
bedtools_closest("-a strand-test-a.bed -b strand-test-b.bed -s")
## break ties
bedtools_closest("-a close-a.bed -b close-b.bed -t first")
## multiple databases
bedtools_closest("-a mq1.bed -b mdb1.bed,mdb2.bed,mdb3.bed -names a,b,c")
## ignoring upstream
bedtools_closest("-a d.bed -b d_iu.bed -D ref -iu")
}
