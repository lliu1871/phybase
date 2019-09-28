\name{write.tree.string}
\alias{write.tree.string}
\title{ Write a tree file }
\description{
  The function writes tree strings to a file in NEXUS or PHYLIP format.
}
\usage{
write.tree.string(X, format = "Nexus", file = "", name = "")
}
\arguments{
  \item{X}{ a vector of tree strings }
  \item{format}{ tree file format }
  \item{file}{ the file name }
  \item{name}{ the species names }
}
\details{
  If name is provided, the function will use name as the species names in the translation block in the NEXUS tree file. Otherwise, the species names will be extracted from the tree strings.}
\value{
  The function returns a tree file in the format of NEXUS or PHYLIP.
}
\author{ Liang Liu \email{lliu@uga.edu} }
\keyword{IO}
\references{
Felsenstein, J. The Newick tree format. \url{http://evolution.genetics.washington.edu/phylip/newicktree.html} 
}
\seealso{ \code{\link{write.subtree}}, \code{\link{read.tree.string}} }
