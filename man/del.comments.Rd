\name{del.comments}
\alias{del.comments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Delete comments }
\description{
  This function deletes comments in the data file.
}
\usage{
del.comments(X)
}
\arguments{
  \item{X}{a vector of strings as the data file is read using scan}
}
\author{ Liang Liu }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304[#0.01],G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
del.comments(treestr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
