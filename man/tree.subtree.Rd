\name{tree.subtree}
\alias{tree.subtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Subtree }
\description{
  The function returns the subtree under the node \code{inode}
}
\usage{
tree.subtree(inode, name, nodematrix)
}
\arguments{
  \item{inode}{ the root node of the subtree }
  \item{name}{ the species names }
  \item{nodematrix}{ the tree node matrix }
}
\value{
  The function returns the tree string of the subtree.}
\author{ Liang Liu \email{lliu@harvard.edu} }
\seealso{ \code{\link{del.node}} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
spname<-read.tree.nodes(treestr)$names
tree.subtree(7,spname,nodematrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
