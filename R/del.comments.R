'del.comments'<-
function(X)
{
      y<-paste(X,collapse="\n")
      if(length(grep("\\[",y))==0) return(X)
      y<-unlist(strsplit(y,split=""))
      left<-which(y=="[")
      right<-which(y=="]")
      if(length(left)!=length(right)) stop("the number of [ != the number of ] in the tree file")
      for(i in length(left):1) y<-y[-(left[i]:right[i])]
      z<-unlist(strsplit(paste(y,sep="",collapse=""),split="\n"))
      z
}
