\name{allfaces}
\alias{allfaces}
\title{All Faces of a Convex Polyhedron}
\usage{
allfaces(hrep)
}
\description{
List all the nonempty faces of a convex polyhedron, giving for each
the dimension, the active set of constraints, and a relative interior point.
}
\arguments{
  \item{hrep}{H-representation of convex polyhedron (see details).}
}
\details{
  See \code{cddlibman.pdf} in the \code{doc} directory of this package,
  especially Sections 1 and 2.

  This function lists all nonempty faces of a convex polyhedron
  given by the H-representation given by the matrix \code{hrep}.
  Let
  \preformatted{
      l <- hrep[ , 1]
      b <- hrep[ , 2]
      v <- hrep[ , - c(1, 2)]
      a <- (- v)
  }

  Then the convex polyhedron in question is the set of
  points \code{x} satisfying
  \preformatted{
      axb <- a \%*\% x - b
      all(axb <= 0)
      all(l * axb == 0)
  }

  A nonempty \emph{face} of a convex polyhedron \eqn{P} is the subset
  of \eqn{P} that is the set of points over which some linear function
  achieves its
  maximum over \eqn{P}.  Note that \eqn{P} is a face of \eqn{P} and appears
  in the list of faces.  By definition the empty set is also a face, but
  is not listed.   These two faces are said to be \emph{improper}, the
  other faces are \emph{proper}.

  A face in the listing is characterized by the set of constraints that
  are \emph{active}, i. e., satisfied with equality, on the face.

  The \emph{relative interior} of a convex set its its interior considered
  as a subset of its affine hull.  The relative interior of
  a one-point set is that point.  The relative interior of a multi-point
  convex set is the union of open line segments \eqn{(x, y)} with endpoints
  \eqn{x} and \eqn{y} in the set.
}
\value{
  a list containing the following components:
  \item{dimension}{list of integers giving the dimensions of the faces.}
  \item{active.set}{list of integer vectors giving for each face the
     set of constraints that are active (satisfied with equality) on
     the face, the integers referring to row numbers of \code{hrep}.}
  \item{relative.interior.point}{list of double or character vectors
     (same type as \code{hrep})
     giving a point in the relative interior of each face.}
}
\section{Rational Arithmetic}{
  The argument \code{hrep} may
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
\seealso{\code{\link{scdd}}, \code{\link{ArithmeticGMP}},
  \code{\link{ConvertGMP}}}
\examples{
hrep <- rbind(c(0, 1,  1,  1, -1),
              c(0, 1,  1, -1, -1),
              c(0, 1, -1, -1, -1),
              c(0, 1, -1,  1, -1),
              c(0, 0,  0,  0,  1))

allfaces(d2q(hrep))
}
\keyword{misc}
