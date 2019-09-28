\name{node.height}
\alias{node.height}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate node height }
\description{
  The function calculates the height of a node. The tree is assumed to be an ultramatric tree. 
}
\usage{
node.height(inode, nodematrix, nspecies)
}
\arguments{
  \item{inode}{ the node number }
  \item{nodematrix}{ the tree node matrix }
  \item{nspecies}{ the number of species in the tree }
}
\value{
The function returns the height of inode.
}
\author{ Liang Liu \email{lliu@uga.edu}}
\keyword{univar}
\examples{
tree.string<-"(((H:4.2,C:4.2):3.1,G:7.3):6.3,O:13.5);"
nodematrix<-read.tree.nodes(tree.string)$nodes
node.height(6,nodematrix,4)

}