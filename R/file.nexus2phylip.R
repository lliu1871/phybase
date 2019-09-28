'file.nexus2phylip'<-
function(inputfolder, outputfolder)
{
	inputfile = list.files(inputfolder)
	outputfile = paste(inputfile,".phy",sep="")
    
	inputfiles<-paste("./",inputfolder,"/",inputfile,sep="")
	outfiles<-paste("./",outputfolder,"/",outputfile,sep="")
    
	for(i in 1:length(inputfiles))
	{
		x<-read.dna.seq(inputfiles[i])
		write.dna.seq(x$seq, name=x$name, file=outfiles[i], format="phylip")
	}
}

