\name{sim.coaltree}
\alias{sim.coaltree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate a coalescence tree }
\description{
This function can simulate a coalescence tree from a single population with parameter theta. The coalescence times in the tree have exponential distributions. \code{theta} is equal to 4uNe where Ne is the effective population size and u is the mutation rate.
}
\usage{
sim.coaltree(nspecies,theta)
}
\arguments{
  \item{nspecies}{ the number of species }
 \item{theta}{ the population parameter }

}
\details{
  theta is the population parameter theta=4N*mu.
}
\value{
  The function returns the simulated coalescence tree.
}
\references{ 
John Wakeley, Coalescent theory: An introduction.
}
\author{ Liang Liu \email{lliu@uga.edu} }
\seealso{  \code{\link{sim.coaltree.sp}} }
\examples{
sim.coaltree(5,theta=0.2)
##[1] "((5:0.55696,(1:0.34858,3:0.34858):0.20838):2.99874,(2:0.97896,4:0.97896):2.57674)"
}
\keyword{ programming }
