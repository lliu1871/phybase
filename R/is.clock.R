'is.clock' <-
function(nodematrix, nspecies,threshold)
{

   clock<-TRUE
   if(dim(nodematrix)[1]==-9)
  	clock<-FALSE
   else{
	for(i in (nspecies+1):(2*nspecies-1)){
		leftson<-nodematrix[i,2]
		rightson<-nodematrix[i,3]
		leftlength<-nodematrix[leftson,4]+ .node.height(leftson,nodematrix,nspecies)
		rightlength<-nodematrix[rightson,4] + .node.height(rightson,nodematrix,nspecies)
		if(abs(leftlength-rightlength)>threshold){
			clock<-FALSE
			break
		}
	}
   }
		
   return(clock)
}

.node.height <- function(inode,nodematrix,nspecies)
{
    if(inode<=nspecies)
    height<-0
    if(inode>nspecies)
    {
        son1<-nodematrix[inode,2]
        height<- .node.height(son1,nodematrix,nspecies)+nodematrix[son1,4]
    }
    return(height)
}

