\name{scdd}
\alias{scdd}
\title{Go between H-representation and V-representation of convex polyhedron}
\usage{
scdd(input, adjacency = FALSE, inputadjacency = FALSE,
    incidence = FALSE, inputincidence = FALSE, roworder = c("lexmin",
    "maxindex", "minindex", "mincutoff", "maxcutoff", "mixcutoff", "lexmax",
    "randomrow"), keepinput = c("maybe", "TRUE", "FALSE"),
    representation = c("H", "V"))
}
\description{
Calculate V-representation (convex hull of points and directions)
of convex polytope given H-representation (intersection of half spaces)
or vice versa.
}
\arguments{
  \item{input}{either H-representation or V-representation of
      convex polyhedron (see details).}
  \item{adjacency}{if \code{TRUE} produce adjacency list of output
      generators.}
  \item{inputadjacency}{if \code{TRUE} produce adjacency list of input
      generators.}
  \item{incidence}{if \code{TRUE} produce incidence list of output
      generators with respect to input generators.}
  \item{inputincidence}{if \code{TRUE} produce incidence list of input
      generators with respect to output generators.}
  \item{roworder}{during the computation, take input rows in the specified
      order.  The default \code{"lexmin"} is usually o. k.  The option
      \code{"maxcutoff"} might be efficient if the input contains many
      redundant inequalities or generators.}
  \item{keepinput}{if \code{"TRUE"} or \code{"maybe"} and an adjacency
      or incidence list involving the input is requested, save the input.}
  \item{representation}{if \code{"H"}, then \code{input} is
      an H-representation, otherwise a V-representation.  May also be
      obtained from a \code{"representation"} attribute of \code{input},
      if present.}
}
\details{
  See \code{cddlibman.pdf} in the \code{doc} directory of this package,
  especially Sections 1 and 2.

  Both representations are (in R) matrices, the first two rows are
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
  \item{output}{An H-representation if input was V-representation
  and vice versa.}
  \item{input}{The argument \code{input}, if requested.}
  \item{adjacency}{The adjacency list, if requested.}
  \item{inputadjacency}{The input adjacency list, if requested.}
  \item{incidence}{The incidence list, if requested.}
  \item{inputincidence}{The input incidence list, if requested.}
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
\seealso{\code{\link{ArithmeticGMP}}, \code{\link{ConvertGMP}},
  \code{\link{validcdd}}, \code{\link{makeH}}}
\examples{
d <- 4
# unit simplex in H-representation
qux <- makeH(- diag(d), rep(0, d), rep(1, d), 1)
print(qux)
# unit simplex in V-representation
out <- scdd(qux)
print(out)
# unit simplex in H-representation
# note: different from original, but equivalent
out <- scdd(out$output)
print(out)

# add equality constraint
quux <- addHeq(1:d, (d + 1) / 2, qux)
print(quux)
out <- scdd(quux)
print(out)

# use some options
out <- scdd(quux, roworder = "maxcutoff", adjacency = TRUE)
print(out)

# convex hull
np <- 50
x <- matrix(rnorm(d * np), ncol = d)
foo <- cbind(0, cbind(1, x))
out <- scdd(d2q(foo), inputincidence = TRUE, representation = "V")
boundies <- sapply(out$inputincidence, length) > 0
sum(boundies) ## number of points on boundary of convex hull
}
\keyword{misc}
