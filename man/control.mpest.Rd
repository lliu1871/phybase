\name{control.mpest}
\alias{control.mpest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ generate a control file for mpest}
\description{
  This function can generate a control file for mpest
}
\usage{
control.mpest(genetreefile, ngene, randomseed=-1, nrun, speciesnames, outputfile)
}
\arguments{
  \item{genetreefile}{the gene tree file}
  \item{ngene}{the number of genes}
  \item{randomseed}{the default is -1; otherwise, a random seed will be generated}
  \item{nrun}{the number of runs; each run has a different starting point, and mp-est will find the tree with the maximum likelihood score across all runs}
  \item{speciesnames}{the names of species}
 \item{outputfile}{the name of the control file}
}
\author{ Liang Liu }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
