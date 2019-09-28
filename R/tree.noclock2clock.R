'tree.noclock2clock' <-
function(inode, treematrix, nspecies)
{
	if(inode > nspecies)
	{
		son1<-treematrix[inode,2]
   		son2<-treematrix[inode,3]
		treematrix <- tree.noclock2clock(son1,treematrix,nspecies)
		treematrix <- tree.noclock2clock(son2,treematrix,nspecies)

        leftheight<-treematrix[son1,4]+node.height(son1,treematrix,nspecies)
		rightheight<-treematrix[son2,4]+node.height(son2,treematrix,nspecies)
 		leftratio<-(leftheight+rightheight)/2/leftheight
		rightratio<-(leftheight+rightheight)/2/rightheight

   		leftsonnodes <- .offspring.nodes(son1,treematrix,nspecies)
   		rightsonnodes<- .offspring.nodes(son2,treematrix,nspecies)

        treematrix[leftsonnodes,4]<-treematrix[leftsonnodes,4]*leftratio
		treematrix[rightsonnodes,4]<-treematrix[rightsonnodes,4]*rightratio
    }
    return(treematrix)
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

