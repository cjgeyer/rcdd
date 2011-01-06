\name{ArithmeticGMP}
\alias{ArithmeticGMP}
\alias{qpq}
\alias{qmq}
\alias{qxq}
\alias{qdq}
\alias{qneg}
\alias{qabs}
\alias{qinv}
\alias{qsum}
\alias{qprod}
\alias{qsign}
\alias{qmatmult}
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
qmatmult(x, y)
qsum(x)
qprod(x)
qsign(x)
qneg(x)
qabs(x)
qinv(x)
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
  \code{qmatmult} is matrix multiplication like \code{\%*\%} for ordinary
  numeric; both arguments must be matrices.
  \code{qneg} is vectorizing negation like unary minus for ordinary numeric.
  \code{qabs} is vectorizing absolute value like \code{abs} for ordinary numeric.
  \code{qinv} is vectorizing inversion like \code{1 / x} for ordinary numeric.
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
qmatmult(matrix(c("1", "2", "3", "4"), 2, 2),
    matrix(c("1/1", "1/2", "1/3", "1/4"), 2, 2))
qneg(seq(-3, 3))
}
\keyword{misc}
