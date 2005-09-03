\name{ConvertGMP}
\alias{ConvertGMP}
\alias{d2q}
\alias{q2d}
\alias{q2q}
\alias{z2q}
\title{Convert Between Real, Integer, and GMP Rational}
\description{
  Converts to and from GMP (GNU multiple precision) rational numbers.
  Any size integers in the numerator and denominator are allowed.
}
\usage{
d2q(x)
q2d(x)
q2q(x)
z2q(numer, denom, canonicalize = TRUE)
}
\arguments{
  \item{x,numer,denom}{objects of type \code{"numeric"}
    or \code{"character"}.}
  \item{canonicalize}{if \code{TRUE} (the default) canonicalize (see below).}
}
\value{
  \code{d2q} converts from real to rational,
  \code{q2d} converts from rational to real,
  \code{q2q} canonicalizes (no common factors in numerator and denominator)
    rationals,
  \code{z2q} converts integer numerator and denominator to rational
    canonicalizing if \code{canonicalize = TRUE} (the default).
}
\seealso{\code{\link{ArithmeticGMP}}}
\examples{
d2q(runif(1))
q2d("-123456789123456789987654321/33")
z2q(44, 11)
}
\keyword{misc}
