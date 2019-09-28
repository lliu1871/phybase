'sptree.steac'<-
function(trees, speciesname, taxaname, species.structure, outgroup, method="nj")
{
    	ntree<-length(trees)
    	nspecies<-length(speciesname)
    	ntax<-length(taxaname)
    	dist <- matrix(0, nrow=ntree, ncol=ntax*ntax) 

    	for(i in 1:ntree)
	{
  		ranktree1 <- read.tree.nodes(trees[i])
		thistreetaxa<-ranktree1$names
		ntaxaofthistree<-length(thistreetaxa)
		ranktree<-ranktree1$nodes
		thistreenode<-rep(-1,ntaxaofthistree)

		#find the missing taxa
		for(j in 1:ntaxaofthistree)
		{
			thistreenode[j]<-which(taxaname == thistreetaxa[j])
			if(length(thistreenode[j])==0)
			{
			print(paste("wrong taxaname",thistreetaxa[j],"in gene",i))
			return(0)
			}
		}

		dist1 <- matrix(0, ntax, ntax)
		for (k in 1:(ntaxaofthistree - 1)) {
        		for (j in (k + 1):ntaxaofthistree) {
            		dist1[thistreenode[k], thistreenode[j]] <- .mrca.2nodes(k, j, ranktree)$dist
        		}
    		} 
   
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
    
	if(method == "upgma"){
		sptree <- tree.upgma(dist,speciesname)$treestr
		sptree <- tree.node2name(sptree,speciesname)
	}
     	if(method == "nj") {
		sptree <- nj(dist)
		sptree <- root(sptree, outgroup = outgroup, resolve.root=TRUE)
		sptree <- write.tree(sptree)
	}  
    	return(sptree)
}

.mrca.2nodes <- function(inode,jnode,nodematrix)
{
    ancestor1<-.ancestor(inode,nodematrix)
    ancestor2<-.ancestor(jnode,nodematrix)
    
    for(i in 1:length(ancestor1)){
        s<-FALSE
        for(j in 1:length(ancestor2)){
            if((ancestor2[j]-ancestor1[i])==0){
                s<-TRUE
                break
            }
        }
        if(s)
        break
    }
    mrca<-ancestor1[i]
    
    totallength<-0
    if(i > 1){
        leftlength<-nodematrix[ancestor1[1:(i-1)],4]
        leftlength<-leftlength[leftlength>0]
        totallength<-sum(leftlength)
    }
    if(j>1)
    {
        rightlength<-nodematrix[ancestor2[1:(j-1)],4]
        rightlength<-rightlength[rightlength>0]
        totallength<-totallength+sum(rightlength)
    }
    z<-list(anc=as.integer,dist=as.double)
    z$anc<-mrca
    z$dist<-totallength
    
    z
}

.ancestor <- function(inode,nodematrix)
{
    if(!is.rootedtree(nodematrix))
    {
        warnings("The tree is not rooted!")
    }
    rootnode <- .rootoftree(nodematrix)
    nnodes<-dim(nodematrix)[1]
    
    if(inode == rootnode)
    ancestor<-rootnode
    else
    {
        ancestor<-rep(0,nnodes)
        ancestor[1]<-inode
        i<-1
        while(ancestor[i] != rootnode)
        {
            ancestor[i+1]<-nodematrix[ancestor[i],1]
            i<-i+1
        }
    }
    return(ancestor[ancestor>0])
}

.rootoftree <- function (nodematrix)
{
    if(sum(nodematrix[,1]<0)>1)
    stop("more than two roots!")
    nodes<-sum((nodematrix[,1]<0)*(1:length(nodematrix[,1])))
    return(nodes)
}


