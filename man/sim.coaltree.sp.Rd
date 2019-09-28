\name{sim.coaltree.sp}
\alias{sim.coaltree.sp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ simulate a gene tree from the species tree }
\description{
  The function simulates a gene tree from the species tree using Rannala and Yang's formula
}
\usage{
sim.coaltree.sp(rootnode, nodematrix, nspecies, seq, name)
}
\arguments{
  \item{rootnode}{ the root node of the species tree }
  \item{nodematrix}{ the tree node matrix of the species tree }
  \item{nspecies}{ the number of species }
  \item{seq}{ a vector of number of sequences in each species }
  \item{name}{species names used in the simulated gene tree. the order of the names must be consistent with that in "nodematrix"}
}
\value{
  \item{gt}{the gene tree generated from the species tree}
  \item{height }{the tree height of the gene tree}
}
\references{ Rannala, B. and Z. Yang. 2003. Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci. Genetics 164: 1645-1656. }
\author{ Liang Liu \email{lliu@uga.edu} }
\seealso{  \code{\link{sim.coaltree}} }
\examples{
tree<-"(((H:0.00402#0.01,C:0.00402#0.01):0.00304#0.01,
G:0.00707#0.01):0.00929#0.01,O:0.01635#0.01)#0.01;"
spname<-species.name(tree)
nodematrix<-read.tree.nodes(tree, spname)$nodes
rootnode<-7
##define the vector seq as [2,2,2,2] which means that there are 2 sequences in each species
seq<-rep(2,4)
str<-sim.coaltree.sp(rootnode,nodematrix,4,seq,name=spname)$gt

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
