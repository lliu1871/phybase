\name{is.clock}
\alias{is.clock}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Is a clock tree or not }
\description{
This function checks the tree to see if the branch lengths satisfy the molecular clock assumption. For each node, the lengths of the left lineage and right lineage are compared. If they are not equal to each other and the difference is greater than \code{threshold}, the function will return FALSE. This function does not perform statistical test for the molecular clock assumption. 
}
\usage{
is.clock(nodematrix, nspecies,threshold)
}
\arguments{
  \item{nodematrix}{ the tree node matrix }
  \item{nspecies}{the number of species}
  \item{threshold}{the critical value for the difference between the length of the left decendant lineage and that of the right decendant lieage of an internode. The difference below the threshold is treated as no difference.}
}
\value{
  The function returns TRUE for a clock tree and FALSE for a non-clock tree.}
\author{ Liang Liu \email{lliu@uga.edu} }
\seealso{\code{\link{is.rootedtree}}}
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00705):0.00929,O:0.01635):0.1,W:0.11635);"
nodematrix<-read.tree.nodes(treestr)$nodes

##if the threshold is set to be large, the tree is a clock tree
is.clock(nodematrix,5,0.0001)
##[1] TRUE

##if the threshold is a small number, the tree is not a clock tree.
is.clock(nodematrix,5,0.00001)
##[1] FALSE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
