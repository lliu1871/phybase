\name{sptree.alleletable}
\alias{sptree.alleletable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create a sequence-species relationship }
\description{
This function can create a matrix to present the sequence-species relationship.
}
\usage{
sptree.alleletable(numsgenenodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numsgenenodes}{ number of sequences for each species }
}
\details{
The matrix created by this function can be used as species.structure.}
\author{ Liang Liu }
\examples{
numsgenenodes<-c(1,1,1,1,1,2,2,1,1,1,1,2,3,2,2,2,1,1,1,2,1,8,2,2,2,1,1,1)
species.structure<-sptree.alleletable(numsgenenodes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
