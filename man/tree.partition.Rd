\name{tree.partition}
\alias{tree.partition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ partition a tree }
\description{
  partition a tree. 
}

\usage{
tree.partition(tree,nspecies)
}
\arguments{
  \item{tree}{ the tree node matrix }
 \item{nspecies}{the number of species}
}
\value{
  The function returns a matrix. Each row represents a particular partition of the tree. The position of "1" in the matrix indicates the presence of the corresponding species in the partition. The last number at each row is the frequency of that partition. This function returns the partition matrix for only one tree. }
\author{ Liang Liu }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
tree.partition(nodematrix,5)
#
#     [,1] [,2] [,3] [,4] [,5] [,6]
#[1,]    1    0    1    0    0    1
#[2,]    1    1    1    0    0    1
#[3,]    1    1    1    1    0    1
#
#The last number of each row is the frequency of the corresponding partition. 
#For example, the frequency of the first partition (1 0 1 0 0) is 1. 
#The first partition includes species 1 and 3 
#as indicated by the position of 1 in the partition. 
#Each row represens a partition and its frequency.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
