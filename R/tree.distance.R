'tree.distance'<-
function(tree1,tree2, method="RF", normalize=TRUE)
{
    tree1 <- read.tree(text=tree1)
    tree2 <- read.tree(text=tree2)
    
    name1 <- tree1$tip.label
    name2 <- tree2$tip.label
    taxaname <- c(name1,name2)
    allspecies <- table(taxaname)
    allspeciesname <- names(allspecies)
    delname <- allspeciesname[which(allspecies==1)]
    if(length(allspeciesname) - length(delname) < 3) return(-1)
    
    x<-drop.tip(tree1,tip=delname)
    y<-drop.tip(tree2,tip=delname)
    
    if(length(x$tip.label)<3)   distance <- -1
    if (method=="RF")    distance <- dist.topo(unroot(x),unroot(y))
    if (method=="SC")    distance <- dist.topo(unroot(x),unroot(y),method="score")
    if (normalize)    distance <- distance/((length(allspeciesname)-length(delname)-3)*2)
    
    z=list(distance=distance, numCommonSp = length(allspeciesname)-length(delname))
    z
}


