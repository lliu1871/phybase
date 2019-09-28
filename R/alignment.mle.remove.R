'alignment.mle.remove'<-
function(path_raxml = "raxmlHPC", seqfiles, contreefile)
{
	reftreefile = paste(seqfiles,".reftree",sep="")
	contree = read.tree(contreefile)
	spnames=contree$tip.label

	for(i in 1:length(seqfiles))
	{
		genename=read.dna.seq(seqfiles[i],format="phylip")$name
		delname=spnames[which(is.na(match(spnames,genename)==1))]
		write.tree(drop.tip(contree,tip=delname),reftreefile[i])
	}

	for(i in 1:length(seqfiles))
	{
		command=paste(path_raxml," -s ",seqfiles[i]," -n ", seqfiles[i]," -mGTRGAMMA -p ",floor(runif(1)*3536219+45497),sep="")
		try(system(command))
	}


	for(i in 1:length(seqfiles))
	{
		x=read.tree(reftreefile[i])
		y=read.tree(paste("RAxML_bestTree.",seqfiles[i],sep=""))
		ntip=length(x$tip.label)
		w=order(x$edge[,2])
		ref_length=(x$edge.length[w])[1:ntip]
		ref_length=ref_length[order(x$tip.label)]
		name = x$tip.label[order(x$tip.label)]

		w=order(y$edge[,2])
		fit_length=(y$edge.length[w])[1:ntip]
		fit_length=fit_length[order(y$tip.label)]
	
		ratio=fit_length[1:ntip]/ref_length[1:ntip]

		if(length(which(ratio>5))>0)
		{
			data=read.dna.seq(seqfiles[i],format="phylip")
			delseq=match(name[which(ratio>5)],data$name)
			write.dna.seq(data$seq[-delseq,],data$name[-delseq],file=paste(seqfiles[i],".final",sep=""),format="phylip")
		}else{
			try(system(paste("cp ",seqfiles[i], " ", paste(seqfiles[i],".final",sep=""))))
		}
	}

}

