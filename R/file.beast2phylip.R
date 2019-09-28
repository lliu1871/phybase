'file.beast2phylip'<-
function(beastfile)
{
	x = scan(beastfile,what="character",sep="\n")

	spname=x[grep("taxon id=",x)]
	spname=gsub("\t\t<taxon id=\"","",spname)
	spname=gsub("\"/>","",spname)
	spname=gsub("\">","",spname)

	gene=grep("nchar",x)

	for(i in 1:length(gene))
	{	
		if(i==length(gene))
		{
		y = x[gene[i]:length(x)]
		y = y[1:grep("/alignment",y)]	
		}else{
		y = x[gene[i]:(gene[i+1]-1)]
		}
		z=unlist(strsplit(y[1],split="ntax="))
		z=unlist(strsplit(z,split="nchar="))
		ntax=0
		nchar=as.numeric(unlist(strsplit(z[3],split=" "))[1])
		w = rep("",length(spname))

		for(j in 1:length(spname))
		{
			if(length(grep(paste("\"",spname[j],"\"",sep=""),y))>0) 
			{
				ntax = ntax + 1
				w[j]=paste(spname[j],y[grep(paste("\"",spname[j],"\"",sep=""),y)+1],collapse="",sep="")
				ww = unlist(strsplit(w[j],split=""))
				ww = ww[(length(ww)-nchar+1):length(ww)]
				ww = toupper(ww)
				if(sum(ww=="A")+sum(ww=="C")+sum(ww=="G")+sum(ww=="T") == 0)
				{
					w[j] = ""
					ntax = ntax - 1
				}
			}
		}
		m = w[w!=""]
		write(c(paste(ntax,nchar),m),paste("gene",i,sep=""))
	}
	return (1)
}
