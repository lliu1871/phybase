'sim.coal.mpest'<-function(mpest_tree, ngenetree)
{
	nodes = read.tree.nodes(mpest_tree)
	spname = nodes$names
	nspecies = length(spname)
	rootnode = dim(nodes$nodes)[1]
	new = tree.noclock2clock(rootnode,nodes$nodes,length(nodes$names))
	new[,5] = new[,4]/nodes$nodes[,4]*2

	genetree=rep("",ngenetree)
	for(i in 1:ngenetree)
	{
		genetree[i]=sim.coaltree.sp(rootnode, new, nspecies, seq=rep(1,nspecies), spname)$gt
	}
	genetree
}

