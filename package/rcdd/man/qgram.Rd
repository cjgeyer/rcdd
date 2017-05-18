\name{qgram}
\alias{qgram}
\title{GMP Rational Gram-Schmidt}
\description{
   Find Orthogonal Basis
}
\usage{
qgram(x, remove.zero.vectors = TRUE)
}
\arguments{
  \item{x}{matrix of type \code{"numeric"}
    or \code{"character"}.  If \code{"numeric"} are converted to rational
    using \code{\link{d2q}}.  Columns are considered vectors in space
    of dimension \code{nrow(x)}.}
  \item{remove.zero.vectors}{logical.}
}
\value{
  an matrix of the same form as \code{x} whose columns are orthogonal
  and span the same vector subspace as the columns of \code{x}.
  Since making the columns unit vectors in the L2 sense could require
  irrational numbers, the columns are made unit vectors in the L1 sense
  unless they are zero vectors (which, of course, cannot be normalized).

  If \code{remove.zero.vectors == TRUE}, then zero vectors are removed
  from the result, so the columns of the result form a basis of the
  subspace.
}
\seealso{\code{\link{ConvertGMP}}}
\examples{
foo <- cbind(c("1", "1", "0", "0", "0"),
             c("2", "1", "0", "0", "0"),
             c("3", "1", "0", "0", "0"),
             c("1", "2", "3", "4", "5"))
qgram(foo)
}
\keyword{misc}
