\name{Subset}
\alias{all.intersect}
\alias{all.union}
\alias{maximal}
\title{subsets and maximal sets}
\usage{
all.intersect(sets, pow2)
all.union(sets, pow2)
maximal(sets, pow2)
}
\description{
Given a list of positive integer vectors representing sets,
return a vector of all pairwise intersections (\code{all.intersect}),
return a vector of all pairwise unions (\code{all.union}),
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
  For \code{all.intersect} or \code{all.union} a
  list of length \code{choose(length(sets), 2)}
  giving all pairwise intersections (resp. unions) of elements of \code{sets}.
  For \code{maximal} a logical vector of the same length as \code{sets}
  indicating the maximal elements.

  Note: \code{all.intersect} and \code{all.union} run over the pairs
  in the same order so they can be matched up.
}
\keyword{misc}
