'tree.brlens'<-
function(tree)
{
    names<-species.name(tree)
    dist<-rep(0,length(names))
    for(j in 1:length(names)) dist[j] <- .dist.tip.root(tree, names[j])
    
    x<-read.tree(text=tree)$edge.length
    list("summary"=summary(x), "sd"=sd(x), "molClock"=sd(dist), "tip.root.dist"=dist)
}

.dist.tip.root<-function(tree,tipname)
{
    nodes<-read.tree.nodes(tree)$node
    name<-read.tree.nodes(tree)$name
    tipnumber<-which(name==tipname)
    root<-dim(nodes)[1]
    length<-0
    father<-tipnumber
    while(father != root)
    {
        length<-length+nodes[father,4]
        father<-nodes[father,1]
    }
    length
}

