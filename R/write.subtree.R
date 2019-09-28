'write.subtree'<-
function (inode, nodematrix, taxaname, root, print.support=FALSE) 
{
    if (nodematrix[inode, 2] > 0) {
	if(nodematrix[inode, 1] == -8)
	{
		if(print.support){
			son1 <- nodematrix[inode, 2]
        		son2 <- nodematrix[inode, 3]
			son3 <- nodematrix[inode, 4]
        		x <- paste("(", write.subtree(son1, nodematrix, taxaname, root, print.support), ":", nodematrix[son1, 4], ",", write.subtree(son2, nodematrix, taxaname, root, print.support), ":", nodematrix[son2, 4], ",", write.subtree(son3, nodematrix, taxaname, root, print.support), ":", nodematrix[son3, 4],")", "[",nodematrix[inode,7],"]",sep = "")
		}else{
			son1 <- nodematrix[inode, 2]
        		son2 <- nodematrix[inode, 3]
			son3 <- nodematrix[inode, 4]
        		x <- paste("(", write.subtree(son1, nodematrix, taxaname, root, print.support), ":", nodematrix[son1, 4], ",", write.subtree(son2, nodematrix, taxaname, root, print.support), ":", nodematrix[son2, 4], ",", write.subtree(son3, nodematrix, taxaname, root, print.support), ":", nodematrix[son3, 4],")", sep = "")
		}

	}else{
		if(print.support){
			son1 <- nodematrix[inode, 2]
        		son2 <- nodematrix[inode, 3] 
        		x <- paste("(", write.subtree(son1, nodematrix, taxaname, 
            root, print.support), ":", nodematrix[son1, 4], ",", write.subtree(son2, 
            nodematrix, taxaname, root, print.support), ":", nodematrix[son2, 
            4], ")", "[",nodematrix[inode,7],"]", sep = "")

		}else{
        		son1 <- nodematrix[inode, 2]
        		son2 <- nodematrix[inode, 3] 
        		x <- paste("(", write.subtree(son1, nodematrix, taxaname, 
            root, print.support), ":", nodematrix[son1, 4], ",", write.subtree(son2, 
            nodematrix, taxaname, root, print.support), ":", nodematrix[son2, 
            4], ")", sep = "")
		}
	}
	
    }
    else x <- taxaname[inode]
    if (inode == root) 
        x <- paste(x, ";", sep = "")
    return(x)
}
