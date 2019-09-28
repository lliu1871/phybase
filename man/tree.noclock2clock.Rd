\name{tree.noclock2clock}
\alias{tree.noclock2clock}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert a non-clocklike tree to a clocklike tree}
\description{
 This function converts a non-clocklike tree to a clocklike tree using an ad-hoc approach described in the paper Liu et al 2007.
}
\usage{
tree.noclock2clock(inode, treematrix, nspecies)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inode}{ root of the tree }
  \item{treematrix}{ tree node matrix }
  \item{nspecies}{ the number of species in the tree}
}
\value{
The function returns the tree node matrix of the clocklike tree.}
\references{ ~put references to the literature/web site here ~ }
\author{ Liang Liu}
\examples{
treestr<-"(((H:1,C:3):2,G:6):2,O:10);"
name<-species.name(treestr)
treenode<-read.tree.nodes(treestr,name)$nodes
tree.noclock2clock(7,treenode,4)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
