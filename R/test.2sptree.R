'test.2sptree'<-function(path_mpest="mpest", sptree1, sptree2, genetreefile, ngenetree, nbootstrap)
{
spname = species.name(sptree1)

tree1_mpest = run.mpest(path_mpest, genetreefile=genetreefile, species=spname, sptree=sptree1, ntree=ngenetree)
tree2_mpest = run.mpest(path_mpest, genetreefile=genetreefile, species=spname, sptree=sptree2, ntree=ngenetree)

if(tree1_mpest$loglike > tree2_mpest$loglike){
	teststat = 2*(tree1_mpest$loglike - tree2_mpest$loglike)
	nulltree = tree2_mpest$tree
}else{
	teststat = 2*(tree2_mpest$loglike - tree1_mpest$loglike)
	nulltree = tree1_mpest$tree
}

nulltree = gsub("0.000000","0.000001",nulltree)

tree1_loglike=tree2_loglike=1:nbootstrap

for(i in 1:nbootstrap)
{
	genetrees = sim.coal.mpest(nulltree,ngenetree=ngenetree)
	write(genetrees,"testboot")
	tree1_loglike[i] = run.mpest(path_mpest, genetreefile="testboot", species=spname, sptree=sptree1, ntree=ngenetree)$loglike
	tree2_loglike[i] = run.mpest(path_mpest, genetreefile="testboot", species=spname, sptree=sptree2, ntree=ngenetree)$loglike
}

z = list(teststat=teststat, value = as.numeric, bootstrap = as.vector, nulltree=nulltree)

if(tree1_mpest$loglike > tree2_mpest$loglike){
	z$bootstrap = 2*(tree1_loglike - tree2_loglike)
}else{
	z$bootstrap = 2*(tree2_loglike - tree1_loglike)
}

z$pvalue = sum(z$bootstrap>=teststat)/nbootstrap

d=density(z$bootstrap)
plot(d,main=NA, xlim=c(min(z$bootstrap),max(z$bootstrap,teststat)),xlab="test statistics", ylab="Density")
abline(v=teststat,lty=3)

z
}
