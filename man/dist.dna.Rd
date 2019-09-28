\name{dist.dna}
\alias{dist.dna}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate pairwise distances among DNA sequences }
\description{
  Calculate pairwise distances among DNA sequences. The sites with missing characters are excluded.
}
\usage{
dist.dna(sequences, nst = 0)
}
\arguments{
  \item{sequences}{ DNA sequences }
  \item{nst}{ substitution model. 0:no model, 1:JC }
}
\details{
  If nst=0, the distance is equal to the proportion of sites with different nucleotides.
}
\value{
  The function returns a distance matrix.
  }
\references{ 
Jukes, TH and Cantor, CR. 1969. Evolution of protein molecules. Pp. 21-123 in H. N. Munro, ed. Mammalian protein metabolism. Academic Press, New York. }

\author{ Liang Liu \email{lliu@uga.edu} }
\seealso{ \code{\link{tree.upgma}} }
\examples{
data(dat.finch)
dist.dna(dat.finch$seq,nst=1)
}
\keyword{ programming }
