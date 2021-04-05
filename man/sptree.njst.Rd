\name{sptree.njst}
\alias{sptree.njst}
\title{calculate the sptree.njst tree}
\description{
This function can estimate species trees from a set of unrooted gene trees 
}
\usage{
sptree.njst(genetrees, spname, taxaname, species.structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genetrees}{ a set of unrooted gene trees}
  \item{spname}{names of species}
  \item{taxaname}{ names of taxa}
  \item{species.structure}{the taxaname-spname table}
}
\author{ Liang Liu }
\examples{
sptree<-"(A:0.4,(B:0.3,(C:0.2,(D:0.1,E:0.1):0.1):0.1):0.1);"

spname<-species.name(sptree)
nspecies<-length(spname)
rootnode<-9
nodematrix<-read.tree.nodes(sptree,spname)$node
seq<-rep(1,nspecies)
species.structure<-matrix(0,nspecies,nspecies)
diag(species.structure)<-1

##population size, theta
nodematrix[,5]<-0.1
ngene<-5
genetree<-rep("",ngene)

##generate gene trees
for(i in 1:ngene)
{
	genetree[i]<-sim.coaltree.sp(rootnode,nodematrix,nspecies,seq,spname)$gt
}
	
##construct the sptree.njst tree
sptree.njst(genetree,spname, spname, species.structure)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
