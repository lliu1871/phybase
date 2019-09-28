\name{test.submodel.valid.raxml}
\alias{test.submodel.valid.raxml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{model validation of the substitution models}
\description{
 The function includes chi-squres goodness-of-fit test for validating the GTRGAMMA model. There are three tests, (1) the chi-square test for base frequencies, (2) the chi-square test for double-nucleotide frequencies, (3) the test for site patterns and the test statistics is sum(|exp_freq-obs_freq|).
}
\usage{
test.submodel.valid.raxml(inputfile, path_raxml, path_seqgen, nbootstrap=100)
}
\arguments{
	\item{inputfile}{the sequence file}
	\item{path_raxml}{the full path of the raxml binary}
	\item{path_seqgen}{the path of the seq-gen binary}
	\item{nbootstrap}{The number of bootstrap replicates in the site pattern test.}
}
\author{ Liang Liu \email{lliu@uga.edu} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
