\name{tree.name2node}
\alias{tree.name2node}
\title{ Replace species names by their node numbers }
\description{
This function replaces the species names in the tree string with their node numbers. 
}
\usage{
tree.name2node(treestr,name="")
}
\arguments{
  \item{treestr}{ the tree string }
  \item{name}{ the species names }
}
\details{
  If species names are not given, the function will use the sorted species names in the tree string. 
}
\value{
  The function returns the tree string with the species names replaced by the node numbers.
}
\author{ Liang Liu \email{lliu@uga.edu} }
\seealso{\code{\link{tree.node2name}} }
\keyword{univar}
\examples{
treestr<-"(((H:4.2,C:4.2):3.1,G:7.3):6.3,O:13.5);"
name<-c("H","G", "C", "O")
tree.name2node(treestr,name)

}