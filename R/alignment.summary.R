'alignment.summary'<-
function(seqfile,format="phylip")
{
	data = read.dna.seq(seqfile,format=format)
	seq = tolower(data$seq)

	nspecies=length(data$name)
	missing_prop=1-(sum(seq=="a")+ sum(seq=="c")+sum(seq=="g")+sum(seq=="t"))/length(seq)
	missing_species=apply(seq,1,function(x) 1-(sum(x=="a")+ sum(x=="c")+sum(x=="g")+sum(x=="t"))/length(x))

	names(missing_species) = data$name
	seq_length=dim(seq)[2]

	for(i in 1:nspecies) if(missing_species[i] == 1) warning(paste("The sequence of species", names(missing_species[i]), "is completely missing"))
	z = list(ntaxa = nspecies, seqlength = seq_length, missing=missing_prop, missing.taxa = missing_species)
	z
}
