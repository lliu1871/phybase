\name{tree.distance}
\alias{tree.distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tree distance}
\description{
  This function calculates the distance between two trees. Two trees are pruned to have the same set of species.
}
\usage{
tree.distance(tree1,tree2, method="RF", normalize=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree1}{a tree string}
  \item{tree2}{a tree string}
  \item{method}{RF: RF distance, SC: branch score distance}
  \item{normalize}{normalized by (2 * the number of internal branches)}
}
\value{
  It calculates the distance of two trees.
}
\author{ Liang Liu }
\examples{
data(dat.unrootedtree)
tree.distance(dat.unrootedtree[1], dat.unrootedtree[2])
}
\keyword{ programming }
