`tree.plot` <-
function(tree)
{
	plot.phylo(read.tree(text=tree))   
}
