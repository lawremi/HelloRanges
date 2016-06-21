\name{argparsing}
\alias{argparsing}
\title{Argument parsing details}

\description{
  HelloRanges uses \CRANpkg{docopt} for parsing the argument string
  passed as the \code{cmd} argument to functions like
  \code{bedtools_intersect}. \command{bedtools} has its own style of
  argument formatting. Here we document the subtle differences.
}

\details{
  Here are the specific differences:
  \itemize{
    \item{
      \CRANpkg{docopt} requires that multi-character arguments are
      prefixed by two hyphens, e.g., \option{--bed}. However,
      \command{bedtools} expects only a single hyphen. It turns out
      \CRANpkg{docopt} is robust to the single-hyphen case, except for
      the first argument. Since the typical convention is to first
      indicate the file, e.g., \option{-a} or \option{-i}, this
      incompatibility does not often arise in practice.
    }
    \item{
      \CRANpkg{docopt} does not allow values of argument to be
      space-separated, while \command{bedtools} often expects space
      separation for multi-valued arguments.  As a compromise,
      HelloRanges expects the values to be comma-separated. Thus,
      \option{-b b.bed c.bed} needs to be \option{-b b.bed,c.bed}.
    }
    \item{
      Most shells support nested commands within parentheses, e.g.,
      \option{-b < (grep foo file.bed)}, but \CRANpkg{docopt} does not
      support that. Instead, nested commands should be enclosed in
      double quotes, e.g., \option{-b < "grep foo file.bed"}. Such
      constructs are supported via \code{\link{pipe}}.
    }
  }
}

\author{Michael Lawrence}