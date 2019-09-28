\name{read.dna.seq}
\alias{read.dna.seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read sequences from files }
\description{
  The function reads sequences from files in the nexus or phylip format.
}
\usage{
read.dna.seq(file="", format="nexus")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{the input file name }
  \item{format}{nexus or phylip}
}
\value{
  \item{seq }{sequences}
  \item{gene}{partitions on the sequences. Each partition represents a gene or a locus.}
}
\author{ Liang Liu }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
