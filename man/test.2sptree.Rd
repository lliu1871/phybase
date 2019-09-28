\name{test.2sptree}
\alias{test.2sptree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{testing if two species trees are significantly different}
\description{
This function is testing if two species trees are significantly different from each other.
}
\usage{
test.2sptree(path_mpest="mpest", sptree1, sptree2, genetreefile, ngenetree, nbootstrap)
}
\arguments{
	\item{path_mpest}{the full path of the mp-est binary}
	\item{sptree1}{the species tree 1}
	\item{sptree2}{the species tree 2}
	\item{genetreefile}{the gene tree file}
	\item{ngenetree}{the number of gene trees}
	\item{nbootstrap}{the number of bootstrap replicates}
}
\author{ Liang Liu \email{lliu@uga.edu} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
