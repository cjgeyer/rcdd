\name{validcdd}
\alias{validcdd}
\title{validate and H-representation or V-representation of convex polyhedron}
\usage{
validcdd(x, representation = c("H", "V"))
}
\description{
Validate an H-representation or V-representation of convex polyhedron,
see \code{\link{scdd}} for description of valid representations.
}
\arguments{
  \item{x}{an H-representation or V-representation to be validated.}
  \item{representation}{if \code{"H"},
      validate \code{x} as an H-representation,
      otherwise as a V-representation.  May also be obtained from
      a \code{"representation"} attribute of \code{x}, if present.}
}
\value{
  always \code{TRUE}.  Fails with error message if not a valid object.
}
\seealso{\code{\link{scdd}}}
\keyword{misc}
