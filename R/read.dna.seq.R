`read.dna.seq` <-
function(file="", format="nexus")
{
    	X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
	X <- gsub("\t"," ",X)
	if(tolower(format) == "nexus")
	{
    		i<-grep("dimension",X,ignore.case=TRUE)
        	string<-gsub(".*ntax","",X[i],ignore.case=TRUE)
		string<-gsub("[=;]"," ",string,ignore.case=TRUE)
		string
    		astring<-unlist(strsplit(string,split=" "))
		astring<-astring[astring!=""]
		ntax <- as.numeric(astring[1])

    		i <- grep("matrix", X, ignore.case = TRUE)
    		j <- grep(";",X)
    		j<-j[j>i][1]
    
    		str<-X[(i+1):(j-1)]
    		str<-gsub("\\[.*\\]","",str)
    		str<-str[str!=""]

    		a<-matrix("",ncol=2,nrow=length(str))

    		for(k in 1:length(str))
		{
			string<-unlist(strsplit(str[k],split=" "))
			string<-string[string!=""]
			a[k,1]<-string[1] 
   			a[k,2]<-paste(string[2:length(string)],collapse="",sep="")
    		}

    		round<-length(str)/ntax
    		seq<-matrix("",nrow=ntax,ncol=2)
    		for(i in 1:ntax)
    		{
			seq[i,1]<-a[i,1]
			seq[i,2]<-paste(a[(0:(round-1))*ntax+i,2],collapse="",sep="")
			seq[i,2]<-gsub(" ","",seq[i,2])
			seq[i,2]<-gsub("\t","",seq[i,2])
    		}		
    		i<-grep("charset",X,ignore.case=TRUE)
    		if(length(i)>0)
		{
    			string<-X[i]
    			string<-gsub(" ","",string)
    			string<-gsub(";","",string)
    			string<-gsub(".*=","",string)
   	 		gene<-as.numeric(unlist(strsplit(string,split="-")))
    			gene<-matrix(gene,ncol=2,byrow=TRUE)
    		}  
		else
    			gene<-"NA"

    		ncha<-nchar(seq[1,2])
    		sequence<-matrix("",nrow=ntax,ncol=ncha)
    		for(i in 1:ntax)
			sequence[i,]<-unlist(strsplit(seq[i,2],split=""))
		taxaname <- seq[,1]
	}
	if(tolower(format) == "phylip")
	{
		string <- unlist(strsplit(X[1],split=" "))
		string <- as.numeric(string[string != ""])
		ntax <- string[1]
		ncha <- string[2]
		taxaname <- rep("",ntax)
		string <- matrix(X[2:length(X)],nrow=ntax)
		sequence <- matrix("",nrow=ntax,ncol=ncha)

		for(i in 1:ntax)
		{
			string1 <- unlist(strsplit(string[i,],split=" "))
			string1 <- string1[string1 != ""]
			taxaname[i] <- string1[1]
			sequence[i,]<-unlist(strsplit(paste(string1[2:length(string1)],collapse="",sep=""),split=""))
		}
		gene <- "NA"
		
	}
	if (tolower(format) == "fasta") {
		index = grep(">",X)
		taxaname = gsub(">","",X[index])
		
		index2 = matrix(0,nrow=length(index),ncol=2)
		index2[,1] = index[1:length(index)]+1
		index2[,2] = c(index[2:length(index)],length(X)+1)-1
		
		seq <- rep("",length(taxaname))
		for(i in 1:dim(index2)[1]){
			seq[i] = paste(X[index2[i,1]:index2[i,2]],collapse="",sep="")
		}
	
		ntax <- length(taxaname)
		ncha <- nchar(seq[1])
		sequence <- matrix("", nrow = ntax, ncol = ncha)
		for (i in 1:ntax) {
			sequence[i,] <- unlist(strsplit(seq[i],split=""))
        	}
        	gene <- "NA"
    	}
	taxaname = gsub("\t","",taxaname)
	taxaname = gsub(" ","",taxaname)
	sequence = gsub("\t","",sequence)
	sequence = gsub(" ","",sequence)

    	z<-list(seq=as.matrix, name="", gene=as.matrix)
    	z$seq<-sequence
    	z$gene<-gene
    	z$name<-taxaname
    	z
}
