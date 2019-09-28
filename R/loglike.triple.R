'loglike.triple'<-
function(sptree,spname,dna)
{
	ntaxa<-dim(dna)[1]
	ntriple<-ntaxa*(ntaxa-1)*(ntaxa-2)/6
	triple<-matrix(0,nrow=ntriple,ncol=5)
	par<-matrix(0,nrow=ntriple,ncol=4)
	nodematrix<-read.tree.nodes(sptree,name=spname)$node

	triple <- .triplenumber(dna)

	for(i in (ntaxa+1):(2*ntaxa-2))
	{
		son1<-nodematrix[i,2]
		sonnode1 <- .offspring.species(son1, nodematrix, ntaxa)
		son2<-nodematrix[i,3]
		sonnode2 <- .offspring.species(son2, nodematrix, ntaxa)
		father <- i

		while(father != (2*ntaxa-1))
		{
			son <- father
			father <- nodematrix[father,1]
			if(nodematrix[father,2] == son)
				sonnode3 <- .offspring.species(nodematrix[father,3],nodematrix,ntaxa)
			else
				sonnode3 <- .offspring.species(nodematrix[father,2],nodematrix,ntaxa)

			for(j in 1:length(sonnode1))
				for(k in 1:length(sonnode2))
					for( w in 1:length(sonnode3))
					{	
						par[sonnode1[j]+sonnode2[k]+sonnode3[w]-5,] <- .triplepara(i,father,nodematrix,ntaxa)
					}
		}
	}

	triplep <- apply(par, 1, .tripleProb)
	loglike <- sum(triple*log(t(triplep)))

	loglike
}

.offspring.species <- function(inode,nodematrix,nspecies)
{
    offspring <- .offspring.nodes(inode,nodematrix,nspecies)
    return(offspring[offspring<=nspecies])
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

.tripleProb <- function (para)
{
    prob <-.C("tripleProb", as.double(para[1]), as.double(para[2]), as.double(para[3]), as.double(para[4]), p0=double(1), p1=double(1), p2=double(1), p3=double(1), p4=double(1), PACKAGE="phybase")
    
    p    <- rep(-1,5)
    p[1] <- prob$p0
    p[2] <- prob$p1
    p[3] <- prob$p2
    p[4] <- prob$p3
    p[5] <- prob$p4
    
    return(p)
}

.triplepara <- function(inode,jnode,nodematrix,nspecies)
{
    par<-rep(0,4)
    height1<-node.height(inode, nodematrix, nspecies)
    height2<-node.height(jnode, nodematrix, nspecies)
    
    if(height1 < height2)
    {
        par[1] <- height2-height1
        par[2] <- height1
        par[3] <- nodematrix[jnode,5]
        par[4] <- nodematrix[inode,5]
    }
    else if(height1 > height2)
    {
        par[1] <- height1-height2
        par[2] <- height2
        par[3] <- nodematrix[inode,5]
        par[4] <- nodematrix[jnode,5]
    }
    else
    {
        warnings("something is wrong in triplepara")
    }
    par
}

.triplenumber <- function(dna)
{
    ntaxa<-dim(dna)[1]
    ntriple<-ntaxa*(ntaxa-1)*(ntaxa-2)/6
    triple<-matrix(0,nrow=ntriple,ncol=5)
    
    for(i in 1:(ntaxa-2))
    for(j in (i+1):(ntaxa-1))
    for(k in (j+1):ntaxa){
        seq<-dna[c(i,j,k),]
        triple[i+j+k-5,]<-site.pattern(seq)[,4]
    }
    triple
}