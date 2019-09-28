\name{file.separateGeneData}
\alias{file.separateGeneData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{separate the concatenated sequences into individual genes}
\description{
  This function separates the concatenated sequences into individual genes.
}
\usage{
file.separateGeneData(nexusfile, missing=c("?","-","N","n"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nexusfile}{sequence data file in nexus format with the character block for individual genes}
  \item{missing}{missing characters}
}
\author{ Liang Liu }
\keyword{ programming }
