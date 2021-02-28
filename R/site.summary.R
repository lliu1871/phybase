'site.summary' <- 
function(sequence, is.cds)
{	
	if(is.cds){
		seqsummary<-rep(-1, 13)
		names(seqsummary) <- c("segregating_site", "informative_site", "ambiguous_site", "seqlength", "ntaxa", "gc_12pos", "gc_12_var","gc_3pos","gc_3_var","basefreq_a", "basefreq_c", "basefreq_g", "basefreq_t")
		seqsummary[1]<-.segsite(sequence)
		seqsummary[2]<-.informsite(sequence)
		seqsummary[3]<-.ambisite(sequence)
		seqsummary[4]<-.seqlength(sequence)
		seqsummary[5]<-dim(sequence)[1]
		seqsummary[6]<-.gc_12pos(sequence)
		seqsummary[7]<-.gc_12pos.variance(sequence)
		seqsummary[8]<-.gc_3pos(sequence)
		seqsummary[9]<-.gc_3pos.variance(sequence)
		seqsummary[10:13]<-.basefreq(sequence)
	}else{
		seqsummary<-rep(-1, 11)
		names(seqsummary) <- c("segregating_site", "informative_site", "ambiguous_site", "seqlength", "ntaxa", "gc", "gc_var","basefreq_a", "basefreq_c", "basefreq_g", "basefreq_t")
		seqsummary[1]<-.segsite(sequence)
		seqsummary[2]<-.informsite(sequence)
		seqsummary[3]<-.ambisite(sequence)
		seqsummary[4]<-.seqlength(sequence)
		seqsummary[5]<-dim(sequence)[1]
		seqsummary[6]<-.gc(sequence)
		seqsummary[7]<-.gc.variance(sequence)
		seqsummary[8:11]<-.basefreq(sequence)
	}
	seqsummary
}

.is.informsite <- function(x)
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

.is.ambisite <- function(x)
{
	x <- tolower(x)
	if(sum(x=="a") + sum(x=="c") + sum(x=="g") +sum(x=="t") == length(x)) y <- 0
	else y<-1
	return(y)	
}

.is.segsite <- function(x)
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

.basefreq <- function(sequence)
{
	seq = tolower(sequence)
	freq <- c(sum(seq == "a"),sum(seq == "c"),sum(seq == "g"),sum(seq == "t"))
	return(freq/sum(freq))
}

.gc <- function(sequence)
{
	seq = tolower(sequence)
	freq <- c(sum(seq == "a"),sum(seq == "c"),sum(seq == "g"),sum(seq == "t"))
	return((freq[2]+freq[3])/sum(freq))
}

.gc.variance <- function(sequence)
{
	seq = tolower(sequence)
	result = apply(seq,1,function(x) (sum(x=="c")+sum(x=="g"))/(sum(x=="c")+sum(x=="g")+sum(x=="a")+sum(x=="t")))
	return(var(result))
}

.gc_3pos <- function(sequence)
{
	if((dim(sequence)[2] %% 3)>0) stop("This is not a gene!")
	seq = tolower(sequence)
	index <- seq(3,dim(seq)[2],by=3)
	seq <- seq[,index]
	freq <- c(sum(seq == "a"),sum(seq == "c"),sum(seq == "g"),sum(seq == "t"))
	return((freq[2]+freq[3])/sum(freq))
}

.gc_3pos.variance <- function(sequence)
{
	if((dim(sequence)[2] %% 3)>0) return(-1)
	index <- seq(3,dim(sequence)[2],by=3)
	sequence1 <- sequence[,index]
	seq = tolower(sequence1)
	result = apply(seq,1,function(x) (sum(x=="c")+sum(x=="g"))/(sum(x=="c")+sum(x=="g")+sum(x=="a")+sum(x=="t")))
	return(var(result))
}

.gc_12pos <- function(sequence)
{
	if((dim(sequence)[2] %% 3)>0) stop("This is not a gene!")
	seq = tolower(sequence)
	index <- seq(3,dim(seq)[2],by=3)
	seq <- seq[,-index]
	freq <- c(sum(seq == "a"),sum(seq == "c"),sum(seq == "g"),sum(seq == "t"))
	return((freq[2]+freq[3])/sum(freq))
}

.gc_12pos.variance <- function(sequence)
{
	if((dim(sequence)[2] %% 3)>0) return(-1)
	index <- seq(3,dim(sequence)[2],by=3)
	sequence1 <- sequence[,-index]
	seq = tolower(sequence1)
	result = apply(seq,1,function(x) (sum(x=="c")+sum(x=="g"))/(sum(x=="c")+sum(x=="g")+sum(x=="a")+sum(x=="t")))
	return(var(result))
}

.informsite <- function(sequence)
{
	sum(apply(sequence,2,.is.informsite))
}

.segsite <- function(sequence)
{
	sum(apply(sequence,2,.is.segsite))
}

.ambisite <- function(sequence)
{
	sum(apply(sequence,2,.is.ambisite))
}

.seqlength <- function(sequence)
{
	dim(sequence)[2]
}

.numtaxa <- function(sequence)
{
	x <- dim(sequence)[2]	
	missing <- 0
	for(i in 1:dim(sequence)[1])
	{
		if(sum(sequence[i,]=="-")+sum(sequence[i,]=="N")+sum(sequence[i,]=="n")+sum(sequence[i,]=="?") == x) missing <- missing + 1
	}
	(dim(sequence)[1] - missing)
}