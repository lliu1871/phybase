\name{write.dna.seq}
\alias{write.dna.seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write sequences to a Nexus file }
\description{
  write sequences to a Nexus file.
}
\usage{
write.dna.seq(sequence, name, file = "", format="nexus", 
program="mrbayes",partition=matrix(0,ncol=2,nrow=1), 
clock=0, popmupr=0, ngen=1000000,nrun=1,nchain=1,samplefreq=100,
taxa=as.vector,burnin=1000,gamma="(3,0.02)", 
outgroup=1,outfile="",append = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{ DNA sequences }
  \item{name}{ taxa names }
  \item{file}{ output file }
  \item{program}{either mrbayes or best.}
  \item{format}{nexus or phylip}
  \item{partition}{ each partition corresponds a gene or a locus. }
  \item{clock}{1:clock, 0:no clock}
  \item{popmupr}{for non-clock species tree model}
  \item{ngen}{number of generations }
  \item{nrun}{number of runs}
  \item{nchain}{number of chains}
  \item{samplefreq}{sampling frequency}
 \item{taxa}{species names if best is defined}
 \item{burnin}{burn in}
 \item{outgroup}{the node number of the outgroup}
 \item{outfile}{output file}
 \item{append}{append or not}
 \item{gamma}{parameters in the inverse gamma distribution as the prior of theta.}
}
\author{ Liang Liu }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
