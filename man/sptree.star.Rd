\name{sptree.star}
\alias{sptree.star}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Build a STAR tree }
\description{
  The function can build a STAR tree from a set of gene trees. Although STAR can handle missing sequences, it requires that 
all possible pairs of species (n choose 2) should appear in at least one gene tree. Otherwise, STAR cannot calculate the pairwise distances
among species. 
}
\usage{
sptree.star(trees, speciesname, taxaname, species.structure,outgroup,method="nj") 
}
\arguments{
  \item{trees}{ the gene tree vector }
  \item{speciesname}{species names}
  \item{taxaname}{taxa names}
  \item{species.structure}{ a matrix defining the species-taxa relationship }
  \item{outgroup}{outgroup}
 \item{method}{UPGMA or NJ}
}
\value{
  The function returns a STAR tree.
}
\author{ Liang Liu \email{lliu@uga.edu} }

\examples{
#create three gene trees
treestr<-rep("",4)
treestr[1]<-"((((H:0.00402,C:0.00402):0.00304,G:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"
treestr[2]<-"((((H:0.00402,G:0.00402):0.00304,C:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"
treestr[3]<-"((((O:0.00402,C:0.00402):0.00304,G:0.00706):0.00929,H:0.01635):0.1,W:0.11635);"
treestr[4]<-"((((H:0.00402,C:0.00402):0.00304,G:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"

speciesname<-species.name(treestr[1])
taxaname<-speciesname
species.structure<-matrix(0,ncol=5,nrow=5)
diag(species.structure)<-1

sptree.star(treestr, speciesname, taxaname, species.structure,outgroup="W",method="nj") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
