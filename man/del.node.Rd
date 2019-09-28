\name{del.node}
\alias{del.node}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Delete a node from the tree }
\description{
  This function deletes a node (and its descendant nodes) from the tree.
}
\usage{
del.node(inode, name, nodematrix)
}
\arguments{
  \item{inode}{ the node to be deleted }
  \item{name}{ the species names }
  \item{nodematrix}{ the tree node matrix  }
}
\details{
  The species names are those defined in the original tree before deleting the node \code{inode}. No need to delete the species name of \code{inode}! If \code{inode} is an internode, the whole subtree below \code{inode} will be deleted.
}
\value{
  \item{nodes }{the tree node matrix after deleting \code{inode}}
  \item{treestr }{the tree string of the tree after deleting \code{inode}. }
}
\author{ Liang Liu }

\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
spname<-read.tree.nodes(treestr)$names
nodematrix<-read.tree.nodes(treestr,spname)$nodes
del.node(6,spname,nodematrix)

##unrooted tree
data(dat.unrootedtree)
name<-paste("S",1:29,sep="")
nodematrix<-read.tree.nodes(dat.unrootedtree[1])$nodes
del.node(6,name,nodematrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
