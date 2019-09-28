'file.phylip2nexus'<-
function(inputfolder, outputfolder)
{
	files   <- paste("./",inputfolder,"/",list.files(inputfolder),sep="")
	outfiles<-paste("./",outputfolder,"/", files,".nex",sep="")
    
	for(i in 1:length(files))
	{
		x<-read.dna.seq(files[i],format="phylip")
		write.dna.seq(x$seq, name=x$name, file=outfiles[i])
	}
}