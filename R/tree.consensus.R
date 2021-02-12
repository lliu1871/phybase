'tree.consensus'<-
function(treefile, outfile, rooted=FALSE, sumtreepath="sumtrees.py")
{
    command <- paste(sumtreepath," -f0 -p ", treefile, " -d0 --replace -F newick -o", outfile, " --no-annotations")
    if(rooted) paste(sumtreepath," -f0 -p ", treefile, " -d0 --replace -F newick -o", outfile, " --no-annotations --rooted")
    try(system(command))
    contree <- read.tree(outfile)
    con = write.tree(contree)
    con
}
