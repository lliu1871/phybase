\name{phybase-package}
\alias{PhyBase}
\docType{package}
\title{Basic functions for Phylogenetic trees}
\description{This package provides functions to read, write, manipulate, simulate, estimate, and summarize phylogenetic trees including species trees which contains not only the topology and branch lengths but also population sizes. The input/output functions can read tree files in which trees are presented in parenthetic format. The trees are read in as a string and then transformed to a matrix which describes the relationship of nodes and branch lengths. The nodes matrix provides an easy access for developers to further manipulate the tree, while the tree string provides interface with other phylogenetic R packages such as "ape". The input/output functions can also be used to change the format of tree files between NEXUS and PHYLIP. Some basic functions have already been established in the package for manipulating trees such as deleting and swapping nodes, rooting and unrooting trees, changing the root of the tree. The package includes functions such as "consensus", "coaltime, "popsize", "treedist" for summarizing phylogenetic trees, calculating the coalescence time, population size, and tree distance. The function maxtree, star.sptree, and steac.sptree are built in the package to esimtate the species tree from multiple gene trees. The packages offers function to simulate DNA sequences from gene trees under substitution models.}
\details{
\tabular{ll}{
Package: \tab PhyBase\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2017-09-26\cr
License: \tab GPL (>=2.0.0)\cr
}

}
\author{
Liang Liu

Maintainer: Liang Liu <lliu@uga.edu>

}
\keyword{ programming }
