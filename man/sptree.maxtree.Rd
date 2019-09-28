\name{sptree.maxtree}
\alias{sptree.maxtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Tree}
\description{
  The function computes the Maximum Tree from multiple gene trees. 
}
\usage{
sptree.maxtree(genetreevector,spname,taxaname,species.structure)
}
\arguments{
  \item{genetreevector}{ a vector of gene trees }
  \item{spname}{the species names}
 \item{taxaname}{the names of taxa}
 \item{species.structure}{the correspondence between species and taxa}
}
\value{
  The function returns the node matrix and tree string of the maximum tree. It also returns the species names.}
\author{ Liang Liu \email{lliu@uga.edu} }
\references{
Liu, L. and D.K. Pearl. Species trees from gene trees: reconstructing Bayesian posterior distributions of a species phylogeny using estimated gene tree distributions. Systematic Biology, 2007, 56:504-514. 

Edwards, S.V., L. Liu., and D.K. Pearl. High resolution species trees without concatenation. PNAS, 2007, 104:5936-5941.

}
\examples{
genetreevector<-c("((((H:0.00302,C:0.00302):0.00304,G:0.00605):0.01029,O:0.01635):0.1,W:0.11635);",
"((((H:0.00402,G:0.00402):0.00304,C:0.00705):0.00929,O:0.01635):0.1,W:0.11635);");
species.structure<-matrix(0,5,5)
diag(species.structure)<-1
name<-species.name(genetreevector[1])
sptree.maxtree(genetreevector,name,name,species.structure)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }

