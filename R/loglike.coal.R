'loglike.coal'<-
function(gtree,sptree,taxaname, spname, species.structure, strict=T)
{
	ntax<-dim(species.structure)[2]
	nspecies<-dim(species.structure)[1]
	ntree<-length(gtree)
	if(!is.matrix(sptree))sptree<-read.tree.nodes(sptree,spname)$nodes
	rootofsptree <- .rootoftree(sptree)
	ab<-matrix(0,nrow=2*nspecies-1,ncol=2)
		
	for(k in 1:ntree)
	{
		genetree<-read.tree.nodes(gtree[k],taxaname)$nodes
		b<-rep(-1,2*nspecies-1)

		coal  <- .getcoaltime(genetree,sptree,ntax,nspecies,species.structure)
		coalt <- matrix(-1,nrow=2*nspecies-1,ncol=3)
		coalt <- .getncoal(rootofsptree,sptree,nspecies,species.structure,coal,coalt)

		node<-which(coalt[,1]>1)
		for(i in 1:length(node))
		{
			height<-sptree[node[i],4]
			time<-rep(height,coalt[node[i],2]+1)
			if(coalt[node[i],2]>0) time[1:coalt[node[i],2]]<-sort(coal[which(coal==node[i]),2])
			time1<-rep(0,length(time))
			if(length(time)>1) time1[2:length(time)]<-time[1:(length(time)-1)]
			time<-time-time1
			n<-coalt[node[i],1]:(coalt[node[i],3])
			b[node[i]]<-sum(n*(n-1)*time)
		}

		ab[,1]<-ab[,1]+coalt[,2]
		ab[,2]<-ab[,2]+b
	}

	ab<-cbind(ab,sptree[,5])
	ab<-ab[which(ab[,2]>-1),]
	if(strict){if(sum(ab[,1]) != (ntax-1)*ntree) stop("the total number of coalescence is not equal to ntax-1")}
	loglike<-sum(ab[,1]*log(2/ab[,3])-(ab[,2]/ab[,3]))

    z=list(coal=coal, loglike=loglike)
	return(z)
}


.getncoal <- function(inode,sptree,nspecies, species.structure,coal,coalt)
{
    
    if(inode <= nspecies)
    {
        coalt[inode,1]<-sum(species.structure[inode,])
        coalt[inode,2]<-sum(coal[,1]==inode)
        coalt[inode,3]<-coalt[inode,1]-coalt[inode,2]
        return(coalt)
        
    }
    else
    {
        coalt <- .getncoal(sptree[inode,2],sptree,nspecies,species.structure,coal,coalt)
        coalt <- .getncoal(sptree[inode,3],sptree,nspecies,species.structure,coal,coalt)
        
        coalt[inode,1]<-coalt[sptree[inode,2],3]+coalt[sptree[inode,3],3]
        coalt[inode,2]<-sum(coal[,1]==inode)
        coalt[inode,3]<-coalt[inode,1]-coalt[inode,2]			
        return(coalt)
    }
}

.getcoaltime <- function(genetree,sptree,ntax,nspecies,species.structure)
{
    geneanc<-rep(0,dim(genetree)[1])
    index<-1
    coal<-matrix(-1,nrow=(dim(genetree)[1]-ntax),ncol=2)
    for(i in 1:nspecies)
    {
        genenode<-which(species.structure[i,]==1)
        spancs<-.ancandtime(i,sptree,nspecies)
        
        for(j in 1:length(genenode))
        {
            inodegene<-genenode[j]
            geneanc1<-.ancandtime(inodegene,genetree,ntax)
            
            geneanc2<-rep(0,dim(genetree)[1])
            geneanc2[geneanc1[,1]]<-1
            newanc<-1:sum((geneanc2-geneanc)>0)
            
            coaltime<-geneanc1[newanc,2]
            for(k in 1:length(coaltime))
            {
                if(k != 1){
                    a<-which(spancs[,2]<=coaltime[k])
                    coal[index,1]<-spancs[length(a),1]
                    coal[index,2]<-coaltime[k]-spancs[a[length(a)],2]
                    index<-index+1}
            }
            geneanc[geneanc1[newanc,1]]<-1
        }
        
    }
    return(coal)
}

.ancandtime<-function (inode, nodematrix,nspecies)
{
    if (!is.rootedtree(nodematrix)) {
        warnings("The tree is not rooted!")
    }
    rootnode <- .rootoftree(nodematrix)
    
    nnodes <- dim(nodematrix)[1]
    if (inode == rootnode)
    ancestor <- rootnode
    else {
        ancestor <- matrix(0, nrow=nnodes,ncol=2)
        ancestor[1,1] <- inode
        ancestor[1,2] <- node.height(inode,nodematrix,nspecies)
        i <- 1
        while (ancestor[i] != rootnode) {
            ancestor[i+1,1] <- nodematrix[ancestor[i,1], 1]
            ancestor[i+1,2] <- ancestor[i,2]+nodematrix[ancestor[i,1],4]
            i <- i + 1
        }
    }
    anc<-ancestor[which(ancestor[,1]>0),]
    colnames(anc)<-c("anc","time")
    anc
}

.rootoftree <- function (nodematrix)
{
    if(sum(nodematrix[,1]<0)>1)
    stop("more than two roots!")
    nodes<-sum((nodematrix[,1]<0)*(1:length(nodematrix[,1])))
    return(nodes)
}
