'sptree.star'<-function (trees, speciesname, taxaname, species.structure, outgroup, method = "nj") 
{
    ntree <- length(trees)
    nspecies <- length(speciesname)
    ntax <- length(taxaname)
    dist <- matrix(0, nrow = ntree, ncol = nspecies * nspecies)
    for (i in 1:ntree) {
        ranktree1 <- read.tree.nodes(trees[i])
        thistreetaxa <- ranktree1$names
        ntaxaofthistree <- length(thistreetaxa)
        ranktree <- ranktree1$nodes
        thistreenode <- rep(-1, ntaxaofthistree)
        for (j in 1:ntaxaofthistree) {
            thistreenode[j] <- which(taxaname == thistreetaxa[j])
            if (length(thistreenode[j]) == 0) {
                print(paste("wrong taxaname", thistreetaxa[j], 
                  "in gene", i))
                return(0)
            }
        }
        if (!is.rootedtree(ranktree)) {
            return("gene trees must be rooted trees!")
        }
        a <- rep(0, 2 * ntaxaofthistree - 1)
        ranknode <- .rank.nodes(ranktree, .rootoftree(ranktree),
            ntaxaofthistree, ntax, a)
        dist1 <- matrix(0, ntax, ntax)
        for (j in (ntaxaofthistree + 1):(2 * ntaxaofthistree - 
            1)) {
            son1 <- .offspring.nodes(ranktree[j, 2], ranktree,
                ntaxaofthistree)
            son1 <- son1[son1 <= ntaxaofthistree]
            son2 <- .offspring.nodes(ranktree[j, 3], ranktree,
                ntaxaofthistree)
            son2 <- son2[son2 <= ntaxaofthistree]
            for (k in 1:length(son1)) for (l in 1:length(son2)) {
                dist1[thistreenode[son1[k]], thistreenode[son2[l]]] <- ranknode[j] * 
                  2
            }
        }
	dist1<-dist1+t(dist1)
	dist1[dist1 == 0] <- NA
	dist1<-dist.species(dist1, species.structure)
        dist[i, ] <- as.numeric(dist1)
    }

    dist2 <- matrix(apply(dist, 2, mean, na.rm = TRUE), nspecies, nspecies)
    diag(dist2) <- 0

    if (sum(is.nan(dist2)) > 0) {
        return("missing species for all genes!")
    }
	
    dist <- dist2
    colnames(dist)<-speciesname
    
    if (method == "upgma") {
            sptree <- tree.upgma(dist, speciesname)$treestr
	    sptree <- tree.node2name(sptree,speciesname)
    }
    if (method == "nj") {
            sptree <- nj(dist)
	    sptree <- root(sptree, outgroup=outgroup, resolve.root = TRUE)
	    sptree <- write.tree(sptree)
    }
   
    return(sptree)
}

.rootoftree <- function (nodematrix)
{
    if(sum(nodematrix[,1]<0)>1)
    stop("more than two roots!")
    nodes<-sum((nodematrix[,1]<0)*(1:length(nodematrix[,1])))
    return(nodes)
}

.rank.nodes <- function(treenode, inode, ntaxa, start, rank)
{
    if(inode > ntaxa){
        left<-treenode[inode,2]
        right<-treenode[inode,3]
        rank[inode]<-start
        rank[left]<-start-1
        rank[right]<-start-1
        rank <- .rank.nodes(treenode,left,ntaxa,rank[left],rank)
        rank <- .rank.nodes(treenode,right,ntaxa,rank[right],rank)
        return(rank)
    }else{
        return(rank)
    } 
}

.offspring.nodes <- function(inode,nodematrix,nspecies)
{
    string <- .offspring.nodes.string(inode,nodematrix,nspecies)
    offspringnodes<-as.numeric(unlist(strsplit(string,split=" ")))
    if(nodematrix[inode,1] == -8)
    offspringnodes<-dim(nodematrix)[1]:1
    return(as.numeric(unlist(strsplit(string,split=" "))))
}

.offspring.nodes.string <- function(inode,nodematrix,nspecies)
{
    if(inode<1 | inode>dim(nodematrix)[1]){
        a<-paste("The node number should be between 1 and",dim(nodematrix)[1])
        stop(a)
    }
    if(inode<=nspecies) return(paste(inode))
    if(inode>nspecies){
        son1 <- nodematrix[inode,2]
        son2 <- nodematrix[inode,3]
        str1 <- .offspring.nodes.string(son1,nodematrix,nspecies)
        str2 <- .offspring.nodes.string(son2,nodematrix,nspecies)
        str  <- paste(inode)
        str  <- paste(str,str1)
        return(paste(str,str2))
    }
}



