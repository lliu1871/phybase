\name{sim.dna}
\alias{sim.dna}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate DNA sequences from substitution models }
\description{
  Simulate DNA sequences from a tree using substitution model 
}
\usage{
sim.dna(nodematrix,seqlength,model,kappa=2,rate=c(1,1,1,1,1,1),
frequency=c(1/4,1/4,1/4,1/4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  	\item{nodematrix}{ the tree node matrix }
 	\item{seqlength}{ sequence length }
	\item{model}{1 JC, 2 H2P, 3 HKY, 4 GTR}
 	\item{kappa}{the transition/transversion ratio}
	\item{rate}{the six rates used in GTR model}
	\item{frequency}{frequencies of four types of nucleotides}
}
\value{
  The function returns DNA sequences simulated from the gene tree \code{nodematrix}. The sequences are coded as 1:A, 2:C, 3:G, 4:T.}
\references{ 
Jukes, TH and Cantor, CR. 1969. Evolution of protein molecules. Pp. 21-123 in H. N. Munro, ed. Mammalian protein metabolism. Academic Press, New York. }
\author{ Liang Liu \email{lliu@uga.edu} }
\seealso{ \code{\link{sim.coaltree}} }
\examples{
tree<-"(((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635);"
nodematrix<-read.tree.nodes(tree)$nodes
sim.dna(nodematrix,100, model=2, kappa=4)
}
\keyword{ programming }
