'del.brlens' <-
function(tree)
{
    tree<-gsub("e-","",tree)
    tree<-gsub("[[:digit:]]","",tree)
    tree<-gsub("\\.","",tree)
    tree<-gsub(":","",tree)
    tree 
}

