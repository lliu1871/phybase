'write.seq.phylip'<-
function(sequence, name, length, outfile = "",append=FALSE) 
{
    ntax  <- length(sequence)
    nchar <- nchar(sequence[1])
    
    add <- nchar%%length
    nround <- (nchar - add)/length
         
    string <- paste(ntax, nchar)
    write(string, file = outfile, append=append)
    for(i in 1:nround)
    {
	if(i == 1)
	{
		string <- paste(name,substr(sequence,(i-1)*length+1,i*length),sep="     ")
		write(string,file=outfile, append=T)
	       write("\n", file=outfile, append=T)
	}else
	{
		write(substr(sequence,(i-1)*length+1,i*length), file=outfile, append=T)
	       write("\n", file=outfile, append=T)
	}
    }
    if(add > 0) write(substr(sequence,nround*length+1,nround*length+add), file = outfile, append=T)
}

