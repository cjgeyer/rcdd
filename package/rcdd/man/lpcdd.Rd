\name{lpcdd}
\alias{lpcdd}
\title{linear programming with exact arithmetic}
\usage{
lpcdd(hrep, objgrd, objcon = as(0, typeof(hrep)), minimize = TRUE,
    solver = c("DualSimplex", "CrissCross"))
}
\description{
Solve linear program or explain why it has no solution.
}
\arguments{
  \item{hrep}{H-representation of convex polyhedron (see details) over
      which an affine function is maximized or minimized.}
  \item{objgrd}{gradient vector of affine function.}
  \item{objcon}{constant term of affine function.}
  \item{minimize}{minimize if \code{TRUE}, otherwise maximize.}
  \item{solver}{type of solver.  Use the default unless you know better.}
}
\details{
  See \code{cddlibman.pdf} in the \code{doc} directory of this package,
  especially Sections 1 and 2 and the documentation of the function
  \code{dd_LPSolve} in Section 4.2.

  This function minimizes or maximizes an affine function \code{x}
  maps to \code{sum(objgrd * x) + objcon} over a convex polyhedron
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
}
\value{
  a list containing some of the following components:
  \item{solution.type}{character string describing the solution type.
    \code{"Optimal"} indicates the optimum is achieved.
    \code{"Inconsistent"} indicates the feasible region is empty
    (no points satisfy the constraints, the polyhedron specified
    by \code{hrep} is empty).
    \code{"DualInconsistent"} or \code{"StrucDualInconsistent"} indicates
    the feasible region is unbounded
    and the objective function is unbounded below when \code{minimize = TRUE}
    or above when \code{minimize = FALSE}.}
  \item{primal.solution}{Returned only when \code{solution.type = "Optimal"},
     the solution to the stated (primal) problem.}
  \item{dual.solution}{Returned only when \code{solution.type = "Optimal"},
     the solution to the dual problem, Lagrange multipliers for the primal
     problem.}
  \item{dual.direction}{Returned only when
     \code{solution.type = "Inconsistent"},
     coefficients of the positive combination of original inequalities
     that proves the inconsistency.}
  \item{primal.direction}{Returned only when
     \code{solution.type = "DualInconsistent"}
     or \code{solution.type = "StrucDualInconsistent"},
     coefficients of the linear combination of columns
     that proves the dual inconsistency, also an unbounded direction
     for the primal LP.}
}
\section{Rational Arithmetic}{
  The arguments \code{hrep}, \code{objgrd}, and \code{objcon} may
  have type \code{"character"} in which case their elements are interpreted
  as unlimited precision rational numbers.  They consist of an optional
  minus sign, a string of digits of any length (the numerator),
  a slash, and another string of digits of any length (the denominator).
  The denominator must be positive.  If the denominator is one, the
  slash and the denominator may be omitted.  The \code{cdd} package
  provides several functions (see \link{ConvertGMP} and \link{ArithmeticGMP})
  for conversion back and forth between R floating point numbers and rationals
  and for arithmetic on GMP rationals.
}
\seealso{\code{\link{scdd}}, \code{\link{ArithmeticGMP}},
  \code{\link{ConvertGMP}}}
\examples{
# first two rows are inequalities, second two equalities
hrep <- rbind(c("0", "0", "1", "1", "0", "0"),
              c("0", "0", "0", "2", "0", "0"),
              c("1", "3", "0", "-1", "0", "0"),
              c("1", "9/2", "0", "0", "-1", "-1"))
a <- c("2", "3/5", "0", "0")
lpcdd(hrep, a)

# primal inconsistent problem
hrep <- rbind(c("0", "0", "1", "0"),
              c("0", "0", "0", "1"),
              c("0", "-2", "-1", "-1"))
a <- c("1", "1")
lpcdd(hrep, a)

# dual inconsistent problem
hrep <- rbind(c("0", "0", "1", "0"),
              c("0", "0", "0", "1"))
a <- c("1", "1")
lpcdd(hrep, a, minimize = FALSE)
}
\keyword{misc}
