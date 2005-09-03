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
  \item{a1}{numerical matrix for inequality constraints.  If vector, treated
      as matrix with one row.}
  \item{b1}{numerical right hand side vector for inequality constraints.}
  \item{a2}{numerical matrix for equality constraints.  If vector, treated
      as matrix with one row.}
  \item{b2}{numerical right hand side vector for equality constraints.}
  \item{x}{if not \code{NULL}, a valid H-representation.}
  \item{a}{numerical matrix for constraints.  If vector, treated as matrix
    with one row.  Constraints are equality in \code{addHeq} and inequality
    in \code{addHin}.}
  \item{b}{numerical right hand side vector for constraints.}

  Arguments \code{a1}, \code{b1}, \code{a2}, and \code{b2} may be missing,
  but must be missing in pairs.
  Rows in \code{x}, if any, are added to new rows corresponding to
  the constraints given by the other arguments.
}
\value{
  a cdd object that can be handed to \code{\link{scdd}}.
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
