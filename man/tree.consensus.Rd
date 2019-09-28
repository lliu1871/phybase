\name{tree.consensus}
\alias{tree.consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{consensus tree}
\description{
  The function builds consensus tree. 
}
\usage{
tree.consensus(treefile, outfile, rooted=FALSE, sumtreepath="sumtrees.py")
}
\arguments{
  \item{treefile}{tree file}
  \item{outfile}{the consensus tree file}
  \item{rooted}{rooted or unrooted trees}
  \item{sumtreepath}{the full path of sumtrees.py}
}

\value{
  \item{contree}{the consensus tree with bootstrap support values}
  \item{bsvalue}{bootstrap values}
}
\author{ Liang Liu \email{lliu@uga.edu} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{programming}
