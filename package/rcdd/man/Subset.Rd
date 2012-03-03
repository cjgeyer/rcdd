\name{Subset}
\alias{all.intersect}
\alias{all.union}
\alias{allIntersect}
\alias{allUnion}
\alias{maximal}
\title{subsets and maximal sets}
\usage{
allIntersect(sets, pow2)
allUnion(sets, pow2)
maximal(sets, pow2)
}
\description{
Given a list of positive integer vectors representing sets,
return a vector of all pairwise intersections (\code{allIntersect}),
return a vector of all pairwise unions (\code{allUnion}),
or a vector indicating the sets that are maximal
in the sense of not being a subset of any other set in the list
(\code{maximal}).  If the list contains duplicate sets,
at most one of each class of duplicates is declared maximal.
}
\arguments{
  \item{sets}{a list of vectors of \code{storage.mode} \code{"integer"}.
  (Unlike most R functions we do not coerce real to integer.)}
  \item{pow2}{use hash table of size \code{2^pow2}.  May be missing.}
}
\value{
  For \code{allIntersect} or \code{allUnion} a
  list of length \code{choose(length(sets), 2)}
  giving all pairwise intersections (resp. unions) of elements of \code{sets}.
  For \code{maximal} a logical vector of the same length as \code{sets}
  indicating the maximal elements.

  Note: \code{allIntersect} and \code{allUnion} run over the pairs
  in the same order so they can be matched up.
}
\section{Note}{The functions \code{allIntersect} and \code{allUnion}
were called \code{all.intersect} and \code{all.union} in previous
versions of this package.  The names were changed because the
\code{\link{all}} function was made generic and these function are
not methods of that one.  These functions were originally intended to
be used to find the faces of a convex set using the output of
\code{\link{scdd}} but now the \code{\link{allfaces}} function does a better
job and does it much more efficiently.  Hence these functions have no known
use, but have not been deleted for reasons of backwards compatibility.
}
\keyword{misc}
