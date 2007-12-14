\name{ArithmeticGMP}
\alias{ArithmeticGMP}
\alias{qpq}
\alias{qmq}
\alias{qxq}
\alias{qdq}
\alias{qsum}
\alias{qprod}
\alias{qsign}
\title{GMP Rational Arithmetic}
\description{
  Add, subtract, multiply, or divide one object to/from/by another using
  GMP (GNU multiple precision) rational arithmetic.
  Any size integers in the numerator and denominator are allowed.
}
\usage{
qpq(x, y)
qmq(x, y)
qxq(x, y)
qdq(x, y)
qsum(x)
qprod(x)
qsign(x)
}
\arguments{
  \item{x,y}{objects of type \code{"numeric"}
    or \code{"character"}.  If \code{"numeric"} are converted to rational
    using \code{\link{d2q}}.  Objects must have the same length.}
}
\details{
  \code{qpq} is \dQuote{plus},
  \code{qmq} is \dQuote{minus},
  \code{qxq} is \dQuote{times},
  \code{qdq} is \dQuote{divide}.
  Divide by zero is an error.  There are no rational NA, NaN, Inf.
  \code{qsum} is vectorizing summation like \code{sum} for ordinary numeric.
  \code{qprod} is vectorizing product like \code{prod} for ordinary numeric.
  \code{qsign} is vectorizing sign like \code{sign} for ordinary numeric.
}
\value{
  an object of the same form as \code{x} that is the sum, difference,
  product, quotient, or sign or (for \code{qsum} and \code{qprod})
  a scalar that is the sum or product.
}
\seealso{\code{\link{ConvertGMP}}}
\examples{
qmq("1/3", "1/2")
# note inexactness of floating point representations
qmq("1/5", 1/5)
qdq("1/5", 1/5)
qsum(c("1", "1/2", "1/4", "1/8"))
qprod(c("1", "1/2", "1/4", "1/8"))
qsign(c("-1", "1/2", "1/-4", "1/8"))
}
\keyword{misc}
