\name{del.brlens}
\alias{del.brlens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Delete branch lengths from trees }
\description{
  This function deletes branch lengths from trees.
}
\usage{
del.brlens(tree)
}
\arguments{
  \item{tree}{trees in the newick format }
}
\author{ Liang Liu }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
del.brlens(treestr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
