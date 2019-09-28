'test.hybrid'<-function(path_mpest="mpest", genetreefile, tree1, tree2, nbootstrap=100)
{
	species = tree1$tip.label
	nspecies = length(species)
	trees = read.tree(genetreefile)

	dist=matrix(0,nrow=length(trees),ncol=2)
	for(i in 1:length(trees))
	{
		dist[i,1]=dist.topo(unroot(tree1),unroot(trees[[i]]))
		dist[i,2]=dist.topo(unroot(tree2),unroot(trees[[i]]))
	}
	index=which(dist[,2]<dist[,1])

	teststat = length(index)/length(trees)

	write.tree(trees[-index],"nuc")

	x = run.mpest(path_mpest, "nuc",species, length(trees)-length(index))

	log = x$loglike
	sptree = x$tree

	simtree = sim.coal.mpest(sptree, ngenetree=10000)
	simtree = read.tree(text=simtree)

	dist=matrix(0,nrow=10000,ncol=2)
	for(i in 1:10000)
	{
		dist[i,1]=dist.topo(unroot(tree1),unroot(simtree[[i]]))
		dist[i,2]=dist.topo(unroot(tree2),unroot(simtree[[i]]))
	}

	testboot = 1:nbootstrap
	for(i in 1:nbootstrap)
	{
		d = dist[sample(1:10000,length(trees),replace=TRUE),]
		index=which(d[,2]<d[,1])
		testboot[i] = length(index)/length(trees)
	}
	
	z = list(teststat=as.numeric, pvalue=as.numeric)
	z$teststat = teststat
	z$pvalue = sum(testboot>teststat)/nbootstrap
	z

}
