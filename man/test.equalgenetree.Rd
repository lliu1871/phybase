\name{test.equalgenetree}
\alias{test.equalgenetree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{testing if gene trees are identical}
\description{
 The function is the likelihood ratio test for the concatenation assumption that all gene trees have the same topology.
}
\usage{
test.equalgenetree(path_raxml, inputfolder, nbootstrap)
}
\arguments{
	\item{path_raxml}{the full path of the raxml binary}
	\item{inputfolder}{the folder that contains all gene data. Each gene is a separate file}
	\item{nbootstrap}{The number of bootstrap replicates}
}
\author{ Liang Liu \email{lliu@uga.edu} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
