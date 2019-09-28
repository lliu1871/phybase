\name{run.mpest}
\alias{run.mpest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{run mpest}
\description{
  The function runs mpest 
}
\usage{
run.mpest(path_mpest="mpest", genetreefile, species, sptree="", ntree)
}
\arguments{
	\item{path_mpest}{the path of mpest}
	\item{genetreefile}{the path of the input gene tree file}
	\item{species}{species names}
	\item{sptree}{the species tree. If the species tree is provided, mp-est will fit the branch lengths for the given tree}
	\item{ntree}{number of gene trees}
}
\author{ Liang Liu \email{lliu@uga.edu} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
