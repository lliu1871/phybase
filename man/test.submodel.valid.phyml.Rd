\name{test.submodel.valid.phyml}
\alias{test.submodel.valid.phyml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{model validation of the substitution models}
\description{
 The function includes chi-squres goodness-of-fit test for validating the substitution model. There are three tests, (1) the chi-square test for base frequencies, (2) the chi-square test for double-nucleotide frequencies, (3) the test for site patterns and the test statistics is sum(|exp_freq-obs_freq|).
}
\usage{
test.submodel.valid.phyml(inputfile, path_phyml, path_seqgen, model, gamma=TRUE, pinv=FALSE, nbootstrap=100)
}
\arguments{
	\item{inputfile}{the sequence file}
	\item{path_phyml}{the full path of the phyml binary}
	\item{path_seqgen}{the path of the seq-gen binary}
	\item{model}{the substitution model}
	\item{gamma}{gamma parameter}
	\item{pinv}{proportion of invariant sites}
	\item{nbootstrap}{The number of bootstrap replicates in the site pattern test.}
}
\author{ Liang Liu \email{lliu@uga.edu} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
