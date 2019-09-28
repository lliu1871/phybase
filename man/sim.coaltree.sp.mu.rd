\name{sim.coaltree.sp.mu}
\alias{sim.coaltree.sp.mu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate a gene tree from the non-clock species tree model }
\description{
  The function generates a random gene tree from the species tree under the non-clock species tree model.
}
\usage{
sim.coaltree.sp.mu(sptree, spname, seq, numgenetree,method="dirichlet",alpha=5.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sptree}{ species tree}
 \item{spname}{species names}
  \item{seq}{ the species-sequences struction, i.e., which sequence belongs to which species }
  \item{numgenetree}{ the number of gene trees to be generated }
 \item{alpha}{the parameter in the gamma distribution. see also \code{mutation_exp}}
 \item{method}{either gamma or dirichlet}
}
\value{
  \item{gt }{the simulated gene tree}
  \item{st }{the node matrix of the species tree}
  \item{seqname}{the names of sequences}
}
\author{ Liang Liu }
\examples{
sptree<-"(((A:0.5,B:0.5):1#0.1,C:1.5):1#0.1,D:2.5)#0.1;"
spname<-c("A","B","C","D")
seq<-c(1,1,1,1) #each species has only one sequence.
sim.coaltree.sp.mu(sptree, spname, seq, numgenetree=1,method="dirichlet",alpha=5.0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
