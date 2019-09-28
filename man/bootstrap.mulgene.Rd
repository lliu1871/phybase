\name{bootstrap.mulgene}
\alias{bootstrap.mulgene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bootstrap sequences from multiple loci }
\description{
  The function bootstraps sequence columns for each locus sampled from the original multilocus data. It consists of two step. First, it bootstraps loci. Then it bootstraps sequences for each locus. 
}
\usage{
bootstrap.mulgene(sequence,gene,name,boot,outfile="")
}
\arguments{
  \item{sequence}{data matrix }
  \item{gene}{location of each locus}
  \item{name}{taxa names of sequences}
  \item{boot}{the number of bootstrap samples}
  \item{outfile}{output file}
}
\details{
 In the sequences matrix, the rows are "Taxa" and the columns are "sites".
}
\value{
  The function generates a data file in phylip format.
}
\author{ Liang Liu \email{lliu@uga.edu} }
\seealso{ \code{\link{bootstrap}}}
\examples{
#construct the DNA sequences of three taxa
seq <- matrix("A",ncol=4,nrow=3)
rownames(seq)<-c("taxa1","taxa2","taxa3")
seq[,2]<-"G"
seq[,3]<-"C"
seq[,4]<-"T" 

name<-rownames(seq) #taxa names of the sequences

#construct two loci. The first two nucleotides represent the first locus, 
#while nucleotide 3 and 4 represent the second locus.
gene<-matrix(0,ncol=2,nrow=2)
gene[1,]<-c(1,2)
gene[2,]<-c(3,4)
gene          
bootstrap.mulgene(seq,gene,name,boot=2,outfile="bootdata.txt")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
