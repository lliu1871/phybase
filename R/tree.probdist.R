'tree.probdist'<-
function(treefile, sumtreepath="sumtrees.py")
{
    command<-paste(sumtreepath, treefile, "--trprobs=post4342p")
    try(system(command))
    x<-scan("post4342p",what="character",sep="\n")
    trees<-x[grep("TREE Tree",x)]
    t<-matrix(unlist(strsplit(trees,split="]")),nrow=2)[2,]
    prob<-matrix(unlist(strsplit(trees,split="]")),nrow=2)[1,]
    y<-unlist(strsplit(prob,split=",probability="))
    y<-y[grep("cumulative",y)]
    y<-unlist(strsplit(y,split=",cumulative"))
    prob<-as.numeric(y[-grep("probability",y)])
    try(system("rm -f post4342p"))
    result<-c(tree=as.character, prob=as.numeric)
    result$tree<-t
    result$prob<-prob
    result
}

