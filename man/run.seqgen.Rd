\name{run.seqgen}
\alias{run.seqgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{run seq-gen}
\description{
  The function runs seq-gen 
}
\usage{
run.seqgen(path_seqgen="./Seq-Gen-1.3.4/seq-gen",nsim=1, seed=123, basefreq=rep(0.25,4), rate=rep(1,6), seqlength=1000, gamma=0, inv=0, treefile, saveformat="phylip", outputfile)
}
\arguments{
	\item{path_seqgen}{the path of seqgen}
	\item{nsim}{The number of simulations. All simulated data are saved to the same outputfile}
	\item{seed}{the random seed}
	\item{basefreq}{base frequencies}
	\item{rate}{six relative rates in the rate matrix}
	\item{seqlength}{the length of simulated sequences}
	\item{gamma}{the gamma parameter}
	\item{inv}{the proportion of invariant sites}
	\item{treefile}{the file of the true tree}
	\item{saveformat}{phlip, nexus, or fasta}
	\item{outputfile}{the simulated sequences are saved to outputfile}
}
\author{ Liang Liu \email{lliu@uga.edu} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
