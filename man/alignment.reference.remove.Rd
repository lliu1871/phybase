\name{alignment.reference.remove}
\alias{alignment.reference.remove}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ alignment.reference.remove removes sequences}
\description{
  This function removes sequences whose branch lengths in the fitted reference tree is 5 times greater than the corresponding branch lengths in the reference tree. The new alignments are saved in the files .removed.
}
\usage{
alignment.reference.remove(path_raxml = "raxmlHPC", seqfiles, nconcatgene)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path_raxml}{the full path of program raxml}
  \item{seqfiles}{the input sequence files}
  \item{nconcatgene}{the number of genes used for building the concatenation tree}
}

\author{ Liang Liu }

\keyword{ programming }
