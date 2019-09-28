\name{tree.upgma}
\alias{tree.upgma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ tree.upgma tree}
\description{
  The function computes the upgma tree from multiple gene trees. 
}
\usage{
tree.upgma(dist, name, method="average")
}
\arguments{
  \item{dist}{ a distance matrix }
  \item{name}{the species names}
  \item{method}{the method for recalculate pairwise distances. two options: averge or min. }
}
\value{
  The function returns a tree node matrix, a tree string and species names.}
\author{ Liang Liu \email{lliu@uga.edu} }
\seealso{\code{\link{sptree.maxtree}}}
\examples{
dist<-matrix(runif(25),5,5)
dist<-(dist+t(dist))/2
diag(dist)<-0
tree.upgma(dist,name=c("H","G","C","O","W"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }

