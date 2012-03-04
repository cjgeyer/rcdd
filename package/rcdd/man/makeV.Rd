\name{makeV}
\alias{makeV}
\alias{addVpoints}
\alias{addVrays}
\alias{addVlines}
\title{make V-representation of convex polyhedron}
\usage{
makeV(points, rays, lines, x = NULL)
addVpoints(points, x)
addVrays(rays, x)
addVlines(lines, x)
}
\description{
Construct V-representation of convex polyhedron.
See \code{\link{scdd}} for description of valid representations.
}
\arguments{
  \item{points}{numerical or character matrix for points.
    If vector, treated as matrix with one row.  Each row is one point.}
  \item{rays}{numerical or character matrix for points.
    If vector, treated as matrix with one row.  Each row represents one ray
      consisting of all nonnegative multiples of the vector which is the row.}
  \item{lines}{numerical or character matrix for points.
    If vector, treated as matrix with one row.  Each row represents one line
      consisting of all scalar multiples of the vector which is the row.}
  \item{x}{if not \code{NULL}, a valid V-representation.}
}
\details{
  In \code{makeV} the arguments \code{points} and \code{rays} and \code{lines}
  may be missing.
}
\section{Rational Arithmetic}{
  The input representation may
  have type \code{"character"} in which case its elements are interpreted
  as unlimited precision rational numbers.  They consist of an optional
  minus sign, a string of digits of any length (the numerator),
  a slash, and another string of digits of any length (the denominator).
  The denominator must be positive.  If the denominator is one, the
  slash and the denominator may be omitted.  This package
  provides several functions (see \link{ConvertGMP} and \link{ArithmeticGMP})
  for conversion back and forth between R floating point numbers and rationals
  and for arithmetic on GMP rationals.

  Arguments may be a mix of numeric and character in which case all are
  converted to GMP rationals (character) and the output is GMP rational.
}
\value{
  a valid V-representation that can be handed to \code{\link{scdd}}.
}
\seealso{\code{\link{scdd}}, \code{\link{validcdd}}}
\examples{
d <- 4
n <- 7
qux <- makeV(points = matrix(rnorm(n * d), ncol = d))
out <- scdd(qux)
out$output
}
\keyword{misc}
