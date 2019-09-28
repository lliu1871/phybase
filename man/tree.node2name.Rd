\name{tree.node2name}
\alias{tree.node2name}
\title{ Replace node numbers by species names in a tree string }
\description{
This function replaces node numbers in a tree string by species names. 
}
\usage{
tree.node2name(treestr,name="")
}
\arguments{
  \item{treestr}{ a tree string }
  \item{name}{ species names }
}
\value{
  The function returns the tree string with the node numbers replaced by the species names.
}
\author{ Liang Liu }
\seealso{\code{\link{tree.name2node}} }
\keyword{univar}
\examples{
treestr<-"(((1:4.2,2:4.2):3.1,3:7.3):6.3,4:13.5);"
name<-c("H","C", "G", "O")
tree.node2name(treestr,name)

}