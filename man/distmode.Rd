\name{distmode}
\alias{distmode}
\alias{distmode,CompressedAtomicList-method}
\alias{distmode,SimpleList-method}
\alias{distmode,vector-method}
\alias{distmode,factor-method}

\title{
  Compute the mode of a distribution
}
\description{
  Computes the mode (and \dQuote{antimode}) of a distribution. It is not
  clear whether this is a generally useful statistic, but
  \command{bedtools} defined it, so we did for completeness.
}
\usage{
distmode(x, anti = FALSE)
}
\arguments{
  \item{x}{
    The vector for which the mode is computed.
  }
  \item{anti}{
    Whether to return the value with the least representation, instead
    of the highest.
  }
}
\details{
  There are methods for List subclasses and ordinary
  vectors/factors. The List methods are useful for aggregation.
}
\value{
  The value that is the most (or least) prevalent in the \code{x}.
}
\author{
  Michael Lawrence
}

\seealso{
  Not to be confused with the data \code{\link{mode}} of a vector.

  \code{\link{bedtools_map}} for an example in the context of aggregation.
}

\examples{
distmode(c(1L, 2L, 1L))
}
