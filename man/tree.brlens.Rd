\name{tree.brlens}
\alias{tree.brlens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summarizing branch lengths}
\description{
  This function summarizes the branch lengths of a tree.
}
\usage{
tree.brlens(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{a tree string}
}
\value{
  \item{summary}{summary of the branch lengths}
  \item{sd}{standard deviation of the branch lengths}
  \item{molClock}{standard deviation of tip-root distances}
}
\author{ Liang Liu }
\examples{
data(dat.unrootedtree)
tree.brlens(dat.unrootedtree[1])
}
\keyword{ programming }
