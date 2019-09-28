"sim.coaltree.sp.mu"<-
function(sptree, spname,seq,numgenetree,method="dirichlet",alpha=5.0)
{
	nodematrix<-read.tree.nodes(sptree,spname)$nodes
	rootnode<-dim(nodematrix)[1]
	nspecies<-(rootnode+1)/2
	ntaxa<-sum(seq)
	
	#generate mutation rates
	#nodematrix<-cbind(nodematrix,rep(-100,rootnode))
	if(tolower(method) == "gamma")
		nodematrix <- .mutation_exp(nodematrix,rootnode,rootnode,nspecies,alpha)
	if(tolower(method) == "dirichlet")
		nodematrix[,6] <- .rdirichlet(1,rep(alpha,dim(nodematrix)[1]))*dim(nodematrix)[1]
	if(tolower(method) == "user")
		nodematrix[,6] <- alpha

	index<-1
	seqname<-rep("",ntaxa)
	for(i in 1:nspecies)
		for(j in 1:seq[i])
		{
			if(seq[i] > 1)
				seqname[index]<-paste(spname[i],"s",j,sep="")
			else
				seqname[index]<-spname[i]
			index<-index+1
		}

	speciesmatrix<-matrix(0,nrow=nspecies,ncol=ntaxa)

	index<-1	
	for(i in 1:length(seq))
	{
		for(j in 1:seq[i])
		{
			speciesmatrix[i,index]<-1
			index<-index+1
		}
	}
	
	spnodedepth<-rep(0,2*nspecies-1)
	for(i in 1:(2*nspecies-1))
	{
		spnodedepth[i]<-node.height(i,nodematrix,nspecies)
	}

	treestr<-rep("",numgenetree)
	for(j in 1:numgenetree)
	{
		str<-sim.coaltree.sp(rootnode,nodematrix,nspecies,seq,name=spname)$gt
		genetree<-read.tree.nodes(str,name=seqname)$nodes
		genenodedepth<-rep(0,2*ntaxa-1)
	
		for(i in 1:(2*ntaxa-1))
		{
			genenodedepth[i]<-node.height(i,genetree,ntaxa)
		}	
		coaltree <- .populationMutation(nodematrix,spnodedepth,genetree,genenodedepth,speciesmatrix)
		treestr[j]<-write.subtree(dim(coaltree)[1],coaltree,seqname,dim(coaltree)[1])
	}
	z <- list(gt=as.character, st=as.matrix,seqname=as.character)
    	z$gt <- treestr
    	z$st <- nodematrix
	z$seqname<-seqname
	return(z)
}

.rdirichlet <- function(n,a)
## pick n random deviates from the Dirichlet function with shape
## parameters a
{
    l<-length(a);
    x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE);
    sm<-x%*%rep(1,l);
    x/as.vector(sm);
}

.mutation_exp<-function(sptree,root,inode,nspecies,alpha)
{
    if(inode == root)
    {
        sptree[root,6]<-1.0
        sptree <- .mutation_exp(sptree,root,sptree[root,2],nspecies,alpha)
        sptree <- .mutation_exp(sptree,root,sptree[root,3],nspecies,alpha)
    }
    else
    {
        sptree[inode,6]<-rgamma(1,alpha,(alpha/sptree[sptree[inode,1],6]))
        if(inode > nspecies)
        {
            sptree <- .mutation_exp(sptree,root,sptree[inode,2],nspecies,alpha)
            sptree <- .mutation_exp(sptree,root,sptree[inode,3],nspecies,alpha)
        }
    }
    return(sptree)
}


.populationMutation<-function (sptree, spnodedepth, genetree, genenodedepth, speciesmatrix)
{  
	index<-1

	nspecies<-(dim(sptree)[1]+1)/2
	ntaxa<-(dim(genetree)[1]+1)/2
	genetreenodes <- rep(-1,2*ntaxa)

	for(i in 1:nspecies)
      	{
		seq<-(1:ntaxa)*(speciesmatrix[i,])
		seq<-seq[seq>0]
        	for(inode in 1:length(seq))
        	{        
        		inodegene <- seq[inode];
        		stop=0;
        		while(inodegene != dim(genetree)[1])
			{
				#check if the node is already taken care of
				for(k in 1:index)
              				if(inodegene == genetreenodes[k]) 
					{
						stop<-1
						break
					}
				if(stop == 1) 	break

				#change the branch length of node p		
				genetree[inodegene,4] <- .ChangeBrlen(sptree, spnodedepth, i, genetree, genenodedepth, inodegene)

				#copy p to genetreenode
				genetreenodes[index] <- inodegene
                 		index<-index+1

				#reset inodegene
				inodegene<-genetree[inodegene,1]
			}
		}
	}
	genetree[,4]<-round(genetree[,4],6)
	return (genetree)

}

.ChangeBrlen<-function(sptree, spnodedepth, spnode, genetree, genenodedepth, genenode)
{
	inode <- .FindSpnodeDownGenenode(sptree,spnodedepth,spnode,genenodedepth,genenode)
	jnode <- .FindSpnodeDownGenenode(sptree,spnodedepth,spnode,genenodedepth,genetree[genenode,1])

	if(inode == jnode)
	{
		length <- (genetree[genenode,4]) * sptree[inode,6]
	}
	else
	{
		father <- sptree[inode,1]
		length <- (spnodedepth[father] - genenodedepth[genenode])*sptree[inode,6]
		while(father != jnode)
		{
			inode <- father;
			father <- sptree[father,1] 
			length <- length + (spnodedepth[father] - spnodedepth[inode])*(sptree[inode,6])
		}
		length <- length + (genenodedepth[genetree[genenode,1]] - spnodedepth[father])*sptree[father,6]
  	}
	return(length)		
}

.FindSpnodeDownGenenode<-function(sptree, spnodedepth, spnode, genenodedepth, genenode)
{
    findnode<-spnode;
    root<-dim(sptree)[1]
    
    depth <- genenodedepth[genenode]
    father <- sptree[spnode,1]
    
    if(genenode > (length(genenodedepth)+1)/2)
    {
        while(spnodedepth[father] <= depth)
        {
            if(father == root)
            {
                findnode <- father
                break
            }
            else
            {
                findnode <- father
                father <- sptree[father,1]
            }
        }
    }    
    return (findnode)
}

