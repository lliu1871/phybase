\name{alignment.mle.remove}
\alias{alignment.mle.remove}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ alignment.mle.remove sequences }
\description{
  This function removes the sequences whose branch lengths in the ML tree is 5 times greater than the branch lengths in the reference tree, i.e., the concatenation tree. The new alignments are saved in the files .final.
}
\usage{
alignment.mle.remove(path_raxml = "raxmlHPC", seqfiles, contreefile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path_raxml}{the full path of program raxml}
  \item{seqfiles}{the input sequence files}
  \item{contreefile}{the concatenation tree file}
}

\author{ Liang Liu }

\keyword{ programming }
