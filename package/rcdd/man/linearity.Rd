\name{linearity}
\alias{linearity}
\title{Find implicit linearities in H-representation and V-representation
   of convex polyhedron}
\usage{
linearity(input, representation = c("H", "V"))
}
\description{
Given V-representation (convex hull of points and directions)
or H-representation (intersection of half spaces) of convex polyhedron
find non-linearity generators that can be made linearity without
changing polyhedron
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
  Interface to the function \code{dd_ImpliedLinearityRows} of the
  \code{cddlib} library,
  see \code{cddlibman.pdf} in the \code{doc} directory of this package,
  especially Sections 1 and 2 and page 9.
  See also \code{\link{scdd}} for a description of the way this package
  codes H-representations and V-representations as R matrices.

  A row of a matrix that is an H-representation or V-representation is
  a linearity row if the first element of that row is 1.  The row is
  an implied linearity row if the first element of that row is 0 but
  if it were 1 the convex polyhedron described would be unchanged.

  The interpretation is as follows.  For an H-representation, the
  linearity (given plus implied) determines the affine hull of the
  polyhedron (the smallest translate of a subspace containing it).
  For a V-representation, the linearity (given plus implied) determines
  the smallest affine set (translate of a subspace) contained in the
  polyhedron.
}
\value{
  a numeric vector, the indices of the implied linearity rows.
  (Note: rows that are linearity rows in the input matrix are not
  contained in this vector.)
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
  (3) occurs roughly one time in ten thousand, and (4) has only occured once
  and only with the \code{\link{redundant}} function.  So the R floating point
  arithmetic version does mostly work, but you cannot trust its results unless
  you can check them independently.
}
\seealso{\code{\link{ArithmeticGMP}}, \code{\link{ConvertGMP}},
  \code{\link{validcdd}}, \code{\link{makeH}}}
\examples{
### calculate affine hull
### determined by given + implied linearity rows
qux <- rbind(c(0, 2, 0, 0, 1),
             c(0, 3, 1, 0, 0),
             c(0, 4, 0, 1, 0),
             c(0, -7, -1, -1, 0))
out <- linearity(qux, representation = "H")
print(out)
qux[out, 1] <- 1
redundant(qux, representation = "H")$output

### calculate minimal nonempty face of polyhedral convex cone
### determined by given + implied linearity rows
qux <- rbind(c(0, 0, 0, 0, 1),
             c(0, 0, 1, 0, 0),
             c(0, 0, 0, 1, 0),
             c(0, 0, -1, -1, 0))
out <- linearity(qux, representation = "V")
print(out)
redundant(qux, representation = "V")$output
}
\keyword{misc}
