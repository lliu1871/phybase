'alignment.site.summary'<-
function(seqfile,format="phylip")
{
	data = read.dna.seq(seqfile,format=format)
	seq = tolower(data$seq)
	site = site.summary(seq)
	site
}
is.informsite <- function(x)
{
	x <- tolower(x)
	m <- 0
	if(length(grep("a",x)) > 1 && length(grep("c",x)) > 1) m <- 1
	if(length(grep("a",x)) > 1 && length(grep("g",x)) > 1) m <- 1
	if(length(grep("a",x)) > 1 && length(grep("t",x)) > 1) m <- 1
	if(length(grep("c",x)) > 1 && length(grep("g",x)) > 1) m <- 1
	if(length(grep("c",x)) > 1 && length(grep("t",x)) > 1) m <- 1
	if(length(grep("g",x)) > 1 && length(grep("t",x)) > 1) m <- 1
	return (m)
}

is.ambisite <- function(x)
{
	x <- tolower(x)
	if(sum(x=="a") + sum(x=="c") + sum(x=="g") +sum(x=="t") == length(x)) y <- 0
	else y<-1
	return(y)	
}

is.segsite <- function(x)
{
	x <- tolower(x)
	na <- sum(x=="a")
	nc <- sum(x=="c")
	ng <- sum(x=="g")
	nt <- sum(x=="t")
	m <- 0
	if(na>0 && nc>0) m<-1
	if(na>0 && ng>0) m<-1
	if(na>0 && nt>0) m<-1
	if(nc>0 && ng>0) m<-1
	if(nc>0 && nt>0) m<-1
	if(ng>0 && nt>0) m<-1
	return(m)
}

basefreq <- function(sequence)
{
	freq <- 1:4
	freq[1] <- sum(tolower(sequence) == "a")
	freq[2] <- sum(tolower(sequence) == "c")
	freq[3] <- sum(tolower(sequence) == "g")
	freq[4] <- sum(tolower(sequence) == "t")
	return(freq/sum(freq))
}

gc <- function(sequence)
{
	freq <- 1:4
	freq[1] <- sum(tolower(sequence) == "a")
	freq[2] <- sum(tolower(sequence) == "c")
	freq[3] <- sum(tolower(sequence) == "g")
	freq[4] <- sum(tolower(sequence) == "t")
	return((freq[2]+freq[3])/sum(freq))
}

gc_3pos <- function(sequence)
{
	if((dim(sequence)[2] %% 3)>0) return(-1)
	index <- seq(3,dim(sequence)[2],by=3)
	sequence1 <- sequence[,index]
	freq <- 1:4
	freq[1] <- sum(tolower(sequence1) == "a")
	freq[2] <- sum(tolower(sequence1) == "c")
	freq[3] <- sum(tolower(sequence1) == "g")
	freq[4] <- sum(tolower(sequence1) == "t")
	return((freq[2]+freq[3])/sum(freq))
}

informsite <- function(sequence)
{
	sum(apply(sequence,2,is.informsite))
}

segsite <- function(sequence)
{
	sum(apply(sequence,2,is.segsite))
}

ambisite <- function(sequence)
{
	sum(apply(sequence,2,is.ambisite))
}

seqlength <- function(sequence)
{
	dim(sequence)[2]
}

numtaxa <- function(sequence)
{
	x <- dim(sequence)[2]	
	missing <- 0
	for(i in 1:dim(sequence)[1])
	{
		if(sum(sequence[i,]=="-")+sum(sequence[i,]=="N")+sum(sequence[i,]=="n")+sum(sequence[i,]=="?") == x) missing <- missing + 1
	}
	(dim(sequence)[1] - missing)
}

site.summary <- function(sequence)
{
	seqsummary<-rep(-1, 11)
	names(seqsummary) <- c("segregating_site", "informative_site", "ambiguous_site", "seqlength", "ntaxa", "gc_content", "gc_3pos","basefreq_a", "basefreq_c", "basefreq_g", "basefreq_t")
	seqsummary[1]<-segsite(sequence)
	seqsummary[2]<-informsite(sequence)
	seqsummary[3]<-ambisite(sequence)
	seqsummary[4]<-seqlength(sequence)
	seqsummary[5]<-dim(sequence)[1]
	seqsummary[6]<-gc(sequence)
	seqsummary[7]<-gc_3pos(sequence)
	seqsummary[8:11]<-basefreq(sequence)
	seqsummary
}
