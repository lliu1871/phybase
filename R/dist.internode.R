'dist.internode' <-
function(tree, taxaname)
{
	ntaxa<-length(taxaname)
	nodematrix<-read.tree.nodes(tree,taxaname)$nodes
	if(is.rootedtree(nodematrix)) nodematrix <- .unroottree(nodematrix)
	dist<-matrix(0, ntaxa,ntaxa)
	for(i in 1:(ntaxa-1))
		for(j in (i+1):ntaxa)
		{
		anc1<-.ancestor(i,nodematrix)
		anc2<-.ancestor(j,nodematrix)
		n<-sum(which(t(matrix(rep(anc1,length(anc2)),ncol=length(anc2)))-anc2==0, arr.ind=TRUE)[1,])-3
		if(n==-1) n<-0
		dist[i,j]<-n
		}
	dist<-dist+t(dist)
	z<-list(dist=as.matrix, taxaname=as.vector)
	z$dist<-dist
	z$taxaname<-taxaname
	z
}

.rootoftree <- function (nodematrix)
{
    if(sum(nodematrix[,1]<0)>1)
    stop("more than two roots!")
    nodes<-sum((nodematrix[,1]<0)*(1:length(nodematrix[,1])))
    return(nodes)
}

.unroottree <- function(nodematrix)
{
    
    nodes<-nodematrix
    root <- .rootoftree(nodematrix)
    newroot<-nodes[root,2]
    nodes[nodes[root,3],4] <- nodes[newroot,4]+nodes[nodes[root,3],4]
    nodes[nodes[root,3],1] <- newroot
    nodes[newroot,4] <- nodes[root,3]
    nodes[newroot,1]<--8
    unrootnodes<-nodes[1:(dim(nodematrix)[1]-1),]
    
    return(unrootnodes)
}

.ancestor <- function(inode,nodematrix)
{
	if(!is.rootedtree(nodematrix)){
		warnings("The tree is not rooted!")
	}
    rootnode <- .rootoftree(nodematrix)
	nnodes<-dim(nodematrix)[1]

    if(inode == rootnode) ancestor<-rootnode
	else{
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

