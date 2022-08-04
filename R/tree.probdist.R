'tree.probdist'<-
function(treefile, sumtreepath="sumtrees.py")
{
    command<-paste(sumtreepath, treefile, "-x=post4342p")
    try(system(command))
    x<-scan("post4342p.topologies.trees",what="character",sep="\n")
    trees<-x[grep("TREE ",x)]
    t<-matrix(unlist(strsplit(trees,split="]")),nrow=2)[2,]
    prob<-as.numeric(gsub("\\[&W","",matrix(t,2)[1,]))
    try(system("rm -f post4342p*"))
    result<-c(tree=as.character, prob=as.numeric)
    result$tree<-matrix(t,2)[2,]
    result$prob<-prob
    result
}


