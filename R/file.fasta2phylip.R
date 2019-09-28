'file.fasta2phylip'<-
function(inputfolder, outputfolder)
{
	fastafiles<-paste("./",inputfolder,"/",list.files(inputfolder),sep="")
	phyfile   <-paste("./",outputfolder,"/",list.files(inputfolder),".phy",sep="")
    
	for(i in 1:length(fastafiles))
	{
		data<-read.dna(fastafiles[i], format="fasta")
		write.dna(data,file=phyfile[i])
	}
}
