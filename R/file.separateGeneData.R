'file.separateGeneData' <-
function(nexusfile, missing=c("?","-","N","n"))
{
    data <- read.dna.seq(nexusfile)
    seq <- data$seq
    name <- data$name
    gene <- data$gene
    ngene <- dim(gene)[1]
    outfile<-paste(nexusfile,".gene",1:ngene,sep="")
    
    for(i in 1:ngene)
    {
        sequence<-seq[,gene[i,1]:gene[i,2]]
        seqlength <- dim(sequence)[2]
        nspecies<-dim(sequence)[1]
        m<-1:nspecies
        for(j in 1:nspecies)
        for(sss in 1:length(missing))
        if(sum(sequence[j,]==missing[sss])==seqlength) m[j]<-0
        index<-m[m>0]
        write.dna.seq(sequence[index,], name=name[index], file=outfile[i], format="phylip")
    }
}
