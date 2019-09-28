\name{write.seq.phylip}
\alias{write.seq.phylip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ write concatenated sequences to a file }
\description{
  This function writes concatenated sequences to a file.
}
\usage{
write.seq.phylip(sequence, name, length, outfile = "",append=FALSE) 
}
\arguments{
  \item{sequence}{concatenated sequences as strings}
  \item{name}{species names}
  \item{length}{the length of sequences per line in the output file}
  \item{outfile}{output file}
  \item{append}{FALSE or TRUE}
}
\author{ Liang Liu }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
