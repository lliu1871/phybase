\name{write.subtree}
\alias{write.subtree}
\title{ Write a sub-tree into a string }
\description{
  write a tree or a sub-tree into a string in parenthetical format
}
\usage{
write.subtree(inode, nodematrix,taxaname,root,print.support=FALSE)
}
\arguments{
  \item{inode}{ the root node of a sub-tree }
  \item{nodematrix}{ a tree node matrix }
  \item{taxaname}{taxa names}
  \item{root}{ the root node of a sub-tree}
  \item{print.support}{print out support values if print.support=TRUE }
}
\details{
  If inode is the root of the tree, the function will write the whole tree into a string in parenthetical format. If inode is not the root node, the function will write the sub-tree into a string. The function works for both rooted trees and unrooted trees.
}
\value{
  The function returns a tree string in parenthetical format
}
\author{ Liang Liu \email{lliu@uga.edu} }
\keyword{IO}
\seealso{ \code{\link{write.tree.string}}, \code{\link{read.tree.nodes}}}
\examples{
data(dat.coaltree)
tree<-read.tree.nodes(dat.coaltree$sptree)
tree$nodes
tree$names
write.subtree(7,tree$nodes,tree$names,7)
}