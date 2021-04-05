'sptree.njst'<-function (genetrees, spname, taxaname, species.structure)
{
    ntree <- length(genetrees)
    ntaxa <- length(taxaname)
    nspecies <- length(spname)
    dist <- matrix(0, nrow = ntree, ncol = nspecies * nspecies)
    for (i in 1:ntree) {
        genetree1 <- read.tree.nodes(genetrees[i])
        thistreetaxa <- genetree1$names
        ntaxaofthistree <- length(thistreetaxa)
        thistreenode <- rep(-1, ntaxaofthistree)
        dist1 <- matrix(0, ntaxa, ntaxa)
        for (j in 1:ntaxaofthistree) {
            thistreenode[j] <- which(taxaname == thistreetaxa[j])
            if (length(thistreenode[j]) == 0) {
                print(paste("wrong taxaname", thistreetaxa[j], 
                  "in gene", i))
                return(0)
            }
        }
        dist1[thistreenode, thistreenode] <- dist.internode (genetrees[i], 
            thistreetaxa)$dist

	dist1[dist1 == 0] <- NA
	dist1<-dist.species(dist1, species.structure)
        dist[i, ] <- as.numeric(dist1)
    }

    dist2 <- matrix(apply(dist, 2, mean, na.rm = TRUE), nspecies, nspecies)
    diag(dist2) <- 0

    if (sum(is.nan(dist2)) > 0) {
        return("missing species for all genes!")
    }

    speciesdistance <- dist2
    tree <- write.tree(nj(speciesdistance))
    tree.node2name(tree, name = spname)
}


.unroottree <- function(nodematrix)
{
    nodes<-nodematrix
    root <- .rootoftree(nodematrix)
    newroot<-nodes[root,2]
    nodes[nodes[root,3],4]<-nodes[newroot,4]+nodes[nodes[root,3],4]
    nodes[nodes[root,3],1]<-newroot
    nodes[newroot,4]<-nodes[root,3]
    nodes[newroot,1]<--8
    unrootnodes<-nodes[1:(dim(nodematrix)[1]-1),]

    return(unrootnodes) 
}

.rootoftree <- function (nodematrix)
{
    if(sum(nodematrix[,1]<0)>1)
    stop("more than two roots!")
    nodes<-sum((nodematrix[,1]<0)*(1:length(nodematrix[,1])))
    return(nodes)
}

