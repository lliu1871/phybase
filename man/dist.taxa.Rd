\name{dist.taxa}
\alias{dist.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate all pairwise distances among taxa in the tree}
\description{
  The function computes all pairwise distances among taxa in the tree. 
}
\usage{
dist.taxa(nodematrix, nspecies)
}
\arguments{
  \item{nodematrix}{ the tree node matrix }
  \item{nspecies}{the number of taxa in the tree}
}
\value{
  The function returns a distance matrix.}
\author{ Liang Liu \email{lliu@uga.edu} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00705):0.00929,O:0.01635):0.1,W:0.11635);"
nodematrix<-read.tree.nodes(treestr)$nodes
dist.taxa(nodematrix,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
