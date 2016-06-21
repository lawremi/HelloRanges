\name{pair}

\alias{pair}
\alias{pair,Vector,Vector-method}
\alias{pair,GenomicRanges,GenomicRanges-method}
\alias{pair,SummarizedExperiment,GenomicRanges-method}
\alias{pair,GAlignments,GenomicRanges-method}

\title{
  Pair up two vectors
}
\description{
  Creates a Pairs from two vectors, optionally via a left outer join.
}
\usage{
pair(x, y, ...)
\S4method{pair}{Vector,Vector}(x, y, by = findMatches(x, y), all.x = FALSE,
                   NA.VALUE = y[NA])
}
\arguments{
  \item{x}{
    The \dQuote{first} vector.
  }
  \item{y}{
    The \dQuote{second} vector.
  }
  \item{by}{
    The \code{\link[S4Vectors]{Hits-class}} object that matches up the
    elements into pairs.
  }
  \item{all.x}{
    If \code{TRUE}, keep every member of \code{x}, even if it has no
    hits. The \dQuote{second} component is filled with the \code{NA.VALUE}.
  }
  \item{NA.VALUE}{
    The value to fill holes in \code{y} when \code{all.x} is
    \code{TRUE}.
  }
  \item{\dots}{
    Arguments for methods.
  }
}
\details{
  This might move to \pkg{S4Vectors} at some point. It is distinct from
  simple \code{\link[S4Vectors]{Pairs}} construction, because it
  performs transformations like a left outer join. More options might be
  added in the future.

  For GRanges and other ranged objects, \code{pair} adds \dQuote{.} to
  the seqlevels, because that is the seqname of the missing GRanges.
}
\value{
  A Pairs object
}
\author{
  Michael Lawrence
}

\seealso{
  \code{\link[S4Vectors]{Pairs-class}}, created by this function.
  \code{\link{bedtools_intersect}}, whose \code{loj} argument motivated
  the creation of this function.
}
