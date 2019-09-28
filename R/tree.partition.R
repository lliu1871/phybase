'tree.partition' <-
function(tree,nspecies)
{
   partition<-matrix(0,(nspecies-2),nspecies)
   for(i in (nspecies+1):(2*nspecies-2)){
       group <- .offspring.species(i,tree,nspecies)
       partition[(i-nspecies),group]<-1
   }
   partition<-cbind(partition,matrix(1,dim(partition)[1],1))
   return(partition)
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
