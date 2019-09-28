'file.concatData' <-
function(inputfiles, confile)
{
    for(i in 1:length(inputfiles))
    {
        if(i==1) x=read.dna.seq(inputfiles[i],format="phylip")$name
        else x=c(x,read.dna.seq(inputfiles[i],format="phylip")$name)
    }
    
    name <- names(table(x))
    data <- .concatData(inputfiles,name)
    
    write.seq.phylip(data$seq, name, length=100, outfile = confile)
}

.concatData <- function(file, spname)
{
	seq <- read.dna.seq(file[1],format="phylip")
	name1 <- seq$name
	index<-match(spname,name1)
	data <- seq$seq[index,]
	sequence<-apply(data,1,function(x) paste(x,collapse="",seq=""))
	sequence <- gsub("NA","?",sequence)
	sequence <- gsub(" ","",sequence)
	genes <- matrix(0,length(file),2)
	genes[1,1] <- 1
	genes[1,2] <- nchar(sequence[1])

	if(length(file)>1)
	{
		for(i in 2:length(file))
		{
		print(i)
		seq <- read.dna.seq(file[i],format="phylip")
		name1 <- seq$name
		index<-match(spname,name1)
		data <- seq$seq[index,]
		sequence1<-apply(data,1,function(x) paste(x,collapse="",seq=""))
		sequence1 <- gsub("NA","?",sequence1)
		sequence1 <- gsub(" ","", sequence1)
		sequence <- paste(sequence,sequence1,sep="")
		genes[i,1] <- genes[i-1,2] + 1
		genes[i,2] <- genes[i,1] + nchar(sequence[1]) - 1
		}
	}
	result <- list(seq=as.character, name=as.character, genes=as.matrix)
	result$seq <- toupper(sequence)
	result$name <- spname
	result$genes <- genes
	return(result)
}

