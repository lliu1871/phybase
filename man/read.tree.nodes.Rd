\name{read.tree.nodes}
\alias{read.tree.nodes}
\title{Read tree nodes}
\description{
  Read a tree string in parenthesic format and output tree nodes, species names and whether the tree is rooted
}
\usage{
read.tree.nodes(str, name = "")
}
\arguments{
  \item{str}{ a tree string in the parenthetical format }
  \item{name}{ species names}
}
\details{ This function reads a tree string into a matrix that describes the relationships among nodes and corresponding branch lengths. Each row in the matrix represents a node. The first n rows contain the information of the nodes at the tips of the tree. The order of the first n nodes is identical to the alphabetic order of the species names given by name. If name is null, the names will be extracted from the tree string and the first n nodes are in the same order as the species names appear in the tree string from left to right.  

The numbers after ":" are branch lengths. The numbers after pound signs are population sizes. The numbers after "%" are mutation rates used in the non-clocklike species tree model.
}
\value{
  \item{nodes}{ nodes is a matrix that describes the relationships among nodes and corresponding branch lengths and population sizes if the tree is a species tree. Each row corresponds a node in the tree. The matrix has 5 columns. The first column is the father of the current node. The following columns are left son, right son, branch length, and population size. The value -9 implies that the information does not exist. The last row is the root of the tree. If the tree is unrooted, the first number of the root node is -8, while it is -9 for a rooted tree.}
  \item{names }{ species names in the same order of the first n nodes.}
  \item{root}{TRUE for a rooted tree, FALSE for an unrooted tree.}

}
\author{ Liang Liu \email{lliu@uga.edu} }
\keyword{IO}
\seealso{ \code{\link{read.tree.string}}, \code{\link{species.name}}  }
\examples{

##read an unrooted tree
data(dat.unrootedtree)
tree<-read.tree.nodes(dat.unrootedtree[1])
tree$nodes
tree$names
tree$root

#read a rooted tree
data(dat.coaltree)
tree<-read.tree.nodes(dat.coaltree$sptree)
tree$nodes
tree$names
tree$root

}
