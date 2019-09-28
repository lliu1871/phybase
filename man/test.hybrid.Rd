\name{test.hybrid}
\alias{test.hybrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{testing if there are two species trees}
\description{
The test is based on the proportion of gene trees supporting the alternative species tree.
}
\usage{
test.hybrid(path_mpest, genetreefile, tree1, tree2, nbootstrap=100)
}
\arguments{
	\item{path_mpest}{the full path of the mpest binary}
	\item{genetreefile}{the input gene tree file}
	\item{tree1}{the null species tree}
	\item{tree2}{the alternative species tree}
	\item{nbootstrap}{the number of bootstrap replicates}
}
\author{ Liang Liu \email{lliu@uga.edu} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
