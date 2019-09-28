'alignment.reference.remove'<-
function(path_raxml="raxmlHPC", seqfiles, nconcatgene)
{
	reftreefile = paste(seqfiles, ".reftree",sep="")

	for(i in 1:nconcatgene)
	{
		if(i==1) x = read.dna.seq(seqfiles[i],format="phylip")$name
		else x = c(x, read.dna.seq(seqfiles[i],format="phylip")$name)
	}
	spnames = sort(names(table(x)))
	file.concatData(seqfiles[1:nconcatgene], "concat.alignments.phy")

	command=paste(path_raxml, " -s concat.alignments.phy -n concat.alignments.phy -mGTRGAMMA -p ",floor(runif(1)*3536219+45487),sep="")
	try(system(command))
	contree = read.tree("RAxML_bestTree.concat.alignments.phy")

	for(i in 1:length(seqfiles))
	{
		genename=read.dna.seq(seqfiles[i],format="phylip")$name
		delname=spnames[which(is.na(match(spnames,genename)==1))]
		write.tree(drop.tip(contree,tip=delname),reftreefile[i])
	}

	for(i in 1:length(seqfiles))
	{
		command=paste(path_raxml, " -s ",seqfiles[i]," -n ",reftreefile[i]," -mGTRGAMMA -p ",floor(runif(1)*3536219+45497)," -f e -t ", reftreefile[i],sep="")
		try(system(command))
	}

	for(i in 1:length(reftreefile))
	{
		x=read.tree(reftreefile[i])
		y=read.tree(paste("RAxML_result.",reftreefile[i],sep=""))
		w=order(x$edge[,2])
		ref_length=x$edge.length[w]
		w=order(y$edge[,2])
		fit_length=y$edge.length[w]
		ntip=length(x$tip.label)
		ratio=fit_length[1:ntip]/ref_length[1:ntip]
		if(length(which(ratio>5))>0)
		{
			data=read.dna.seq(seqfiles[i],format="phylip")
			delseq=match(x$tip.label[which(ratio>5)],data$name)
			write.dna.seq(data$seq[-delseq,],data$name[-delseq], file=paste(seqfiles[i],".removed",sep=""), format="phylip")
		}else{
			try(system(paste("cp", seqfiles[i], paste(seqfiles[i],".removed",sep=""))))	
		}
	}

}

