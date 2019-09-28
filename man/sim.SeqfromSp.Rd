\name{sim.SeqfromSp}
\alias{sim.SeqfromSp}
\title{simulate DNA sequences from a species tree }
\description{
The function simulates sequences from a species tree. The function assumes that seq-gen has been installed
}
\usage{
sim.SeqfromSp(sptree, spname, ntaxasp, ngene, theta=0, noclock=0, 
simsequence=1, murate="Dirichlet",alpha=5, seqlength=100, rate=c(1,1,1,1,1,1), frequency=c(1/4,1/4,1/4,1/4), 
outfile, format="phylip", concat=TRUE)}
\arguments{
  \item{sptree}{A species tree which must be a rooted tree. }
  \item{spname}{species names}
  \item{ntaxasp}{a vector of the number of individuals in each species}
\item{ngene}{number of genes}
\item{theta}{population size}
\item{noclock}{0: clocklike species tree 1: nonclocklike species tree}
\item{simsequence}{1: simulate sequences and gene trees, 0: simulate gene trees} 
\item{murate}{distribution of mutation rates}
\item{alpha}{the shape parameter of dirichlet distribution}
\item{seqlength}{the number of nucleotides along the sequences}
\item{rate}{rates}
\item{frequency}{nucleotide frequency}
\item{outfile}{the full path of the output file}
\item{format}{either "phylip" or "nexus"}
\item{concat}{save the concatenated sequences or save single-gene sequences as different data in the same file}
}
\value{
  The function writes sequences into a file.
}
\author{ Liang Liu \email{lliu@uga.edu} }
\keyword{programming}

