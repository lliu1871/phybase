'del.node'<-
function (inode, name, nodematrix) 
{
    nodes <- nodematrix
    root <- .rootoftree(nodematrix)
    if (inode == root) 
        stop("Delete the whole tree!")
    if (nodes[inode, 1] == root) {
        if (nodes[root, 3] == inode) 
            newroot <- nodes[root, 2]
        if (nodes[root, 2] == inode) 
            newroot <- nodes[root, 3]
    }
    else {
        newroot <- root
        father <- nodes[inode, 1]
        grandfather <- nodes[father, 1]
        if (nodes[father, 3] == inode) 
            son <- nodes[father, 2]
        if (nodes[father, 2] == inode) 
            son <- nodes[father, 3]
        nodes[inode, 4] <- nodes[inode, 4] + nodes[father, 4]
        nodes[son, 1] <- grandfather
        nodes[son, 4] <- nodes[son, 4] + nodes[father, 4]
        if (nodes[grandfather, 3] == father) 
            nodes[grandfather, 3] <- son
        if (nodes[grandfather, 2] == father) 
            nodes[grandfather, 2] <- son
    }
    z <- list(nodes = matrix(0, dim(nodes)[1], 5), treestr = "")
    z$nodes <- nodes
    z$treestr <- write.subtree(newroot, nodes, name, newroot)
    z
}

.swap.nodes <- function(inode, jnode, name, nodematrix)
{
    nodes<-nodematrix
    ##swap father
    father1<-nodes[inode,1]
    father2<-nodes[jnode,1]
    nodes[inode,1]<-father2
    nodes[jnode,1]<-father1
    
    ##swap son
    if(nodes[father1,2] == inode)
    nodes[father1,2]<-jnode
    else
    nodes[father1,3]<-jnode
    if(nodes[father2,2] == jnode)
  		nodes[father2,2]<-inode
        else
        nodes[father2,3]<-inode
        z <- list(nodes = matrix(0, dim(nodes)[1], dim(nodes)[2]), treestr="")
        z$nodes<-nodes
        z$treestr<-tree.subtree(dim(nodes)[1],name,nodes)
        z 
}

.rootoftree <- function (nodematrix)
{
    if(sum(nodematrix[,1]<0)>1)
    stop("more than two roots!")
    nodes<-sum((nodematrix[,1]<0)*(1:length(nodematrix[,1])))
    return(nodes)
}
