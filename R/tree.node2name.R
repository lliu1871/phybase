'tree.node2name' <-
function(treestr, name="")
{
    str<-treestr
    str<-gsub(" ","",str)    
    if(length(name)==0) stop("please provide names")
    
    if(length(grep(":",str))==0)
    {
    	for(i in length(name):1) 
	{
		str<-gsub(paste("\\(",i,",",sep=""),paste("\\(",name[i],",",sep=""),str)
		str<-gsub(paste(",",i,"\\)",sep=""),paste(",",name[i],"\\)",sep=""),str)
		str<-gsub(paste(",",i,",",sep=""),paste(",",name[i],",",sep=""),str)
	}
    }else{
   	old<-paste(1:length(name),":",sep="")
    	new<-paste(name,"$",sep="")
    	for(i in length(name):1) str<-gsub(old[i],new[i],str)
    	str<-gsub("\\$",":",str)
   }
   str
}
