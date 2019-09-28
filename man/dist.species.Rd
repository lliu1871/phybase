\name{dist.species}
\alias{dist.species}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate pairwise distances among species}
\description{
If some species have multiple taxa, the pairwise distance between two species is equal to the average of the distances between all pairs of taxa in the two species. This functions returns the pairwise distances among species (average over all taxa in the species).
}
\usage{
dist.species(dist, species.structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{ the distance matrix of taxa }
  \item{species.structure}{a matrix with rows representing species and columns representing taxa. 1: the species (row) has the taxon at the corresponding column. see the example.  }
}
\value{
This functions returns the distance matrix of species.
}
\author{ Liang Liu }
\seealso{ See Also as \code{\link{dist.taxa}} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00705):0.00929,O:0.01635):0.1,W:0.11635);"
nodematrix<-read.tree.nodes(treestr)$nodes
dist<-dist.taxa(nodematrix,5)
species.structure<-matrix(0,nrow=2,ncol=5) #2 species and 5 taxa
species.structure[1,]<-c(1,1,1,0,0)        #taxa 1,2,3 belong to the first species
species.structure[2,]<-c(0,0,0,1,1)        #taxa 4,5 belong to the second species
dist.species(dist,species.structure)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
