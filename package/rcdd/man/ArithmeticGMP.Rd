\name{ArithmeticGMP}
\alias{ArithmeticGMP}
\alias{qpq}
\alias{qmq}
\alias{qxq}
\alias{qdq}
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
}
\value{
  an object of the same form as \code{x} that is the sum, difference,
  product, or quotient.
}
\seealso{\code{\link{ConvertGMP}}}
\examples{
qmq("1/3", "1/2")
# note inexactness of floating point representations
qmq("1/5", 1/5)
qdq("1/5", 1/5)
}
\keyword{misc}
