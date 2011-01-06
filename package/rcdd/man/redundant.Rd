\name{redundant}
\alias{redundant}
\title{Eliminate redundant rows of H-representation and V-representation}
\usage{
redundant(input, representation = c("H", "V"))
}
\description{
Eliminate redundant rows from H-representation (intersection of half spaces)
or V-representation (convex hull of points and directions) of convex polytope.
}
\arguments{
  \item{input}{either H-representation or V-representation of
      convex polyhedron (see details).}
  \item{representation}{if \code{"H"}, then \code{input} is
      an H-representation, otherwise a V-representation.  May also be
      obtained from a \code{"representation"} attribute of \code{input},
      if present.}
}
\details{
  See \code{cddlibman.pdf} in the \code{doc} directory of this package,
  especially Sections 1 and 2.

  Both representations are (in R) matrices, the first two columns are
  special.  Let \code{foo} be either an H-representation or
  a V-representation and
  \preformatted{
      l <- foo[ , 1]
      b <- foo[ , 2]
      v <- foo[ , - c(1, 2)]
      a <- (- v)
  }

  In the H-representation the convex polyhedron in question is the set of
  points \code{x} satisfying
  \preformatted{
      axb <- a \%*\% x - b
      all(axb <= 0)
      all(l * axb == 0)
  }

  In the V-representation the convex polyhedron in question is the set of
  points \code{x} for which there exists a \code{lambda} such that
  \preformatted{
      x <- t(lambda) \%*\% v
  }
  where \code{lambda} satisfies the constraints
  \preformatted{
      all(lambda * (1 - l) >= 0)
      sum(b * lambda) == max(b)
  }

  An H-representation or V-representation object can be checked for validity
  using the function \code{\link{validcdd}}.
}
\value{
  a list containing some of the following components:
  \item{output}{The input matrix with redundant rows removed.}
  \item{implied.linearity}{For an H-representation, row numbers of inequality
      constraint rows that together imply equality constraints.  For a
      V-representation, row numbers of rays that together imply lines.}
  \item{redundant}{Row numbers of redundant rows.  Note: this is set
      \code{redset} output by the \code{dd_MatrixCanonicalize} function
      in \code{cddlib}.  It apparently does not consider all rows it deletes
      \dQuote{redundant}.  Redundancy can also be determined from the
      following component.}
  \item{new.position}{Integer vector of length \code{nrow(input)}.  Says
      for each input row which output row it becomes or zero to indicate
      redundant.}
}
\section{Rational Arithmetic}{
  The input representation may
  have type \code{"character"} in which case its elements are interpreted
  as unlimited precision rational numbers.  They consist of an optional
  minus sign, a string of digits of any length (the numerator),
  a slash, and another string of digits of any length (the denominator).
  The denominator must be positive.  If the denominator is one, the
  slash and the denominator may be omitted.  The \code{cdd} package
  provides several functions (see \link{ConvertGMP} and \link{ArithmeticGMP})
  for conversion back and forth between R floating point numbers and rationals
  and for arithmetic on GMP rationals.
}
\section{Warning}{
  If you want correct answers, use rational arithmetic.  If you do not,
  this function may (1) produce approximately correct answers, (2) fail with
  an error, (3) give answers that are nowhere near correct with no error or
  warning, or (4) crash R losing all work done to that point.  In large
  simulations (1) is most frequent, (2) occurs roughly one time in a thousand,
  (3) occurs roughly one time in ten thousand, and (4) has only occured once.
  So the R floating point
  arithmetic version does mostly work, but you cannot trust its results unless
  you can check them independently.
}
\seealso{\code{\link{ArithmeticGMP}}, \code{\link{ConvertGMP}},
  \code{\link{validcdd}}, \code{\link{makeH}}}
\examples{
hrep <- rbind(c(0, 0,  1,  1,  0),
              c(0, 0, -1,  0,  0),
              c(0, 0,  0, -1,  0),
              c(0, 0,  0,  0, -1),
              c(0, 0, -1, -1, -1))

redundant(d2q(hrep), representation = "H")

foo <- c(1, 0, -1)
hrep <- cbind(0, 1, rep(foo, each = 9), rep(foo, each = 3), foo)
print(hrep)
redundant(d2q(hrep), representation = "V")
}
\keyword{misc}
