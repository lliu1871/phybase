\name{read.tree.string}
\alias{read.tree.string}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read tree strings from a tree file }
\description{
  This function reads tree strings in Newick format from a tree file. The output of the function is a vector of tree strings that can be converted to a matrix of nodes by the function \code{\link{read.tree.nodes}}.
}
\usage{
read.tree.string(file = "", format="nexus")
}
\arguments{
  \item{file}{ the tree file that contains trees in Newick format.}
  \item{format}{format = "nexus" or format = "phylip"}
}
\details{
The function can read NEXUS and PHYLIP tree files. It works for other types of tree files as long as the trees in the tree files are in Newick format. This function combining with \code{\link{write.tree.string}} can change the tree file format. 
}
\value{
  \item{tree}{a vector of tree strings.}
  \item{names}{species names.}
  \item{root}{TRUE for rooted trees, FALSE for unrooted trees}
}
\author{ Liang Liu \email{lliu@uga.edu} }
\keyword{IO}
\seealso{ \code{\link{write.tree.string}}, \code{\link{read.tree.nodes}}}
\examples{
##read rooted trees in PHYLIP format
cat("(((H:4.2,C:4.2):3.1,G:7.3):6.3,O:13.5);",file = "phylip.tre", sep = "\n")
tree.string<-read.tree.string("phylip.tre",format="phylip")
}
