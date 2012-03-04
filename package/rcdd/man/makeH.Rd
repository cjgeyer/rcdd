\name{makeH}
\alias{makeH}
\alias{addHeq}
\alias{addHin}
\title{make H-representation of convex polyhedron}
\usage{
makeH(a1, b1, a2, b2, x = NULL)
addHeq(a, b, x)
addHin(a, b, x)
}
\description{
Construct H-representation of convex polyhedron,
set of points \code{x} satisfying
\preformatted{
    a1 \%*\% x <= b1
    a2 \%*\% x == b2
}
see \code{\link{scdd}} for description of valid representations.
}
\arguments{
  \item{a1}{numerical or character matrix for inequality constraints.
    If vector, treated as matrix with one row.}
  \item{b1}{numerical or character right hand side vector for inequality
    constraints.}
  \item{a2}{numerical or character matrix for equality constraints.  If
    vector, treated as matrix with one row.}
  \item{b2}{numerical or character right hand side vector for equality
    constraints.}
  \item{x}{if not \code{NULL}, a valid H-representation.}
  \item{a}{numerical or character matrix for constraints.  If vector,
    treated as matrix with one row.  Constraints are equality
    in \code{addHeq} and inequality in \code{addHin}.}
  \item{b}{numerical or character right hand side vector for constraints.}

  Arguments \code{a1}, \code{b1}, \code{a2}, and \code{b2} may be missing,
  but must be missing in pairs.
  Rows in \code{x}, if any, are added to new rows corresponding to
  the constraints given by the other arguments.
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
  a valid H-representation that can be handed to \code{\link{scdd}}.
}
\seealso{\code{\link{scdd}}, \code{\link{validcdd}}}
\examples{
d <- 4
# unit simplex in H-representation
qux <- makeH(- diag(d), rep(0, d), rep(1, d), 1)
print(qux)
# add an inequality constraint
qux <- addHin(c(1, -1, 0, 0), 0, qux)
print(qux)
# drop a constraint
qux <- qux[- 3, ]
print(qux)
}
\keyword{misc}
