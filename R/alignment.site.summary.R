'alignment.site.summary'<-
function(seqfile,format="phylip",is.cds=FALSE)
{
	data = read.dna.seq(seqfile,format=format)
	seq = tolower(data$seq)
	site = site.summary(seq,is.cds)
	site
}
