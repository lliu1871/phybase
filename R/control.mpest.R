'control.mpest'<-
function(genetreefile, ngene, randomseed=-1, nrun, speciesnames, outputfile)
{
	write(genetreefile, outputfile)
	write("0",outputfile, append=TRUE)
	if(randomseed==-1) write('-1',outputfile, append=TRUE)
	else write(floor(runif(1)*38728212+8847564), outputfile, append=TRUE)
	write(nrun,outputfile,append=TRUE)
	write(paste(ngene,length(speciesnames),sep=" "),outputfile, append=TRUE)
	write(paste(speciesnames,1,speciesnames,sep=" "),outputfile, append=TRUE)
	write("0",outputfile, append=TRUE)
}

