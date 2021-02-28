\name{alignment.site.summary}
\alias{alignment.site.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ summary statistics of sequences }
\description{
  This function calculates the number of segregating sites, informative sites, ambiguous sites (i.e., sites containing non-nucleotide characters) in the alignment.
}
\usage{
alignment.site.summary(seqfile,format="phylip",is.cds=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqfile}{the input sequence files; one gene per file}
  \item{format}{format = phylip, fasta, nexus}
  \item{is.cds}{is it cds?}
}

\author{ Liang Liu }

\keyword{ programming }
