\name{bootstrap}
\alias{bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bootstrap sequences }
\description{
  This function can be used to bootstrap sequences.
}
\usage{
bootstrap(sequence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{sequence matrix.}
}
\details{
 In the sequences matrix, the columns are "Taxa" and the rows are "sites". The function will bootstrap the rows.
}
\value{
  the function returns a sequence matrix with sites randomly sampled from the original matrix with replacement.
}
\author{ Liang Liu }
\examples{
#construct the DNA sequences of three taxa
seq <- matrix("A",ncol=4,nrow=3)
rownames(seq)<-c("taxa1","taxa2","taxa3")
seq[,2]<-"G"
seq[,3]<-"C"
seq[,4]<-"T" 
bootstrap_sample = bootstrap(seq)

data(dat.finch)
bootstrap_sample <- bootstrap(dat.finch$seq)
}
\keyword{ programming }
