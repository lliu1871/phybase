'sim.SeqfromSp' <-
function (sptree, spname, ntaxasp, ngene, theta=0, noclock=0, simsequence=1, murate="Dirichlet",alpha=5, seqlength=100, rate=c(1,1,1,1,1,1), frequency=c(1/4,1/4,1/4,1/4), outfile, format="phylip", concat=TRUE)
{
	if(!is.rootedtree(sptree))
		stop("the species tree must be rooted!")

	nspecies<-length(spname)
	ntaxa <- sum(ntaxasp)

	if(length(ntaxasp) != nspecies)
	{
		stop("The number of species in ntaxasp does not match the number of species in sptree!")
	}

	rootnode<-nspecies*2-1
	taxaname<-rep("", ntaxa)

	index <- 1
	
	for(i in 1:nspecies)
	{
		for(j in 1:ntaxasp[i])
		{
			if(ntaxasp[i] > 1)
				taxaname[index]<-paste(spname[i],"s",j,sep="")
			else
				taxaname[index]<-spname[i]
			index <- index + 1
		}
	}

	species.structure<-matrix(0,nspecies, ntaxa)
	index <- 1
	for(i in 1:nspecies)
	{
		for(j in 1:ntaxasp[i])
		{
			species.structure[i,index]<-1
			index <- index + 1
		}
	}
	

	#simulate genetrees
	nodematrix<-read.tree.nodes(sptree,spname)$nodes
	if(theta > 0)
		nodematrix[,5]<-theta
	tree <- rep("", ngene)

	dnaseq <- matrix("", nrow = ntaxa, ncol = seqlength * ngene)
	partition <- matrix(0, nrow = ngene, ncol = 2)

	for(i in 1:ngene)
	{
		if(noclock == 0)
			tree[i]<-sim.coaltree.sp(.rootoftree(nodematrix),nodematrix, nspecies,ntaxasp, name=spname)$gt
		else if(noclock == 1)
		{			
			tree[i]<-sim.coaltree.sp.mu(sptree, spname, ntaxasp, 1, method=murate,alpha=alpha)$gt
		}
		else
			stop ("noclock could be 0 or 1!")
    }
    genetreefile<-paste(outfile,"_genetree.tre",sep="")

    if(simsequence == 0)
    {
	write(tree, outfile)
    }else{
    if(ngene==1) write(tree, genetreefile)
    else write(paste("[",seqlength,"]",tree,sep=""),genetreefile)
    randomseed <- floor(runif(1)*47362822683 + 7362822683)
    if(tolower(format) == "phylip")
    {
        if(concat) commandline<-paste("seq-gen -mGTR -z", randomseed, " -f", paste(frequency,collapse=",")," -r",paste(rate,collapse=",")," -p ", ngene," -l",paste(ngene*seqlength)," <", genetreefile," > ", outfile)
        else commandline<-paste("seq-gen -mGTR -z", randomseed, " -f", paste(frequency,collapse=",")," -r",paste(rate,collapse=",")," -p ", ngene," -l",paste(seqlength)," <", genetreefile," > ", outfile)
        try(system(commandline))
    }

    if(tolower(format) == "nexus")
    {
        for(i in 1:ngene) partition [i, ] <- c(1+(i-1)*seqlength, seqlength*i)
        commandline<-paste("seq-gen -mGTR -z", randomseed, " -f", paste(frequency,collapse=",")," -r",paste(rate,collapse=",")," -p ", ngene," -l",paste(ngene*seqlength)," <", genetreefile," > ", outfile)
        try(system(commandline))
        seqfile<-read.dna.seq(outfile,format="phylip")
		write.dna.seq(seqfile$seq, name=seqfile$name, file=outfile,partition = partition, format=format, taxa = ntaxasp, program="best", append=FALSE)
    }}
}

.rootoftree <- function (nodematrix)
{
    if(sum(nodematrix[,1]<0)>1)
    stop("more than two roots!")
    nodes<-sum((nodematrix[,1]<0)*(1:length(nodematrix[,1])))
    return(nodes)
}
