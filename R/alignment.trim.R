'alignment.trim'<-
function(path_trimal,inputfolder, outputfolder)
{
	inputfile = list.files(inputfolder)
	output = paste(inputfile,".trimmed",sep="")

	for(i in 1:length(output)) 
	{
		command = paste(path_trimal, " -in ", "./", inputfolder,"/", inputfile[i]," -out ", "./", outputfolder,"/", output[i]," -gappyout", sep="")
		try(system(command))
	}
}
