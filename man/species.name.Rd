\name{species.name}
\alias{species.name}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Species names in a tree string }
\description{
  The function can be used to obtain species names from a tree string.
}
\usage{
species.name(str)
}
\arguments{
  \item{str}{ a tree string in the parenthetical format }
}
\details{
  The function returns the species names. If the tree string contains only the node number instead of species names, the function will return the node numbers.
}
\value{
  The function returns the species names.
}
\author{ Liang Liu \email{lliu@uga.edu} }
\keyword{programming}
\seealso{ \code{\link{read.tree.string}} }
\examples{
tree.string<-"(((H:4.2,C:4.2):3.1,G:7.3):6.3,O:13.5);"
species.name(tree.string)


}