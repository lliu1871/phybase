\name{site.pattern}
\alias{site.pattern}
\title{ Site patterns }
\description{
  The function returns site patterns. 
}
\usage{
site.pattern(seq)
}
\arguments{
  \item{seq}{ DNA sequences with rows representing taxa and columns representing sites }
}
\value{
  The function returns a matrix. Each row in the matrix represents a site pattern and the last number at each row is the frequency of the site pattern appeared in the DNA sequences.  
}
\author{ Liang Liu \email{lliu@uga.edu} }

\examples{
seq<- matrix("A",nrow=4,ncol=5)
seq[1,]<-c("A","A","G","C","C")
seq[2,]<-c("A","G","G","C","C")
seq[3,]<-c("T","A","G","C","C")
seq[4,]<-c("A","A","G","T","T")
site.pattern(seq)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
