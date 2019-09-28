'dist.taxa' <-
function(nodematrix,nspecies)
{
    dist<-matrix(0,nspecies,nspecies)
    for(i in 1:(nspecies-1)){
	for(j in (i+1):nspecies){
		dist[i,j] <- .mrca.2nodes(i,j,nodematrix)$dist
        }
    }
    return(dist+t(dist))
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
    rootnode<-rootoftree(nodematrix)
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

