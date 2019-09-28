`dist.dna` <-
function(sequences,nst=0)
{
    sequences = tolower(sequences)
    nsequences<-dim(sequences)[1]

    dist<-matrix(0,nsequences,nsequences)
    for(i in 1:(nsequences-1)){
        for(j in (i+1):nsequences){
            x = which(sequences[i,]=="a" | sequences[i,]=="c" | sequences[i,]=="g" | sequences[i,]=="t")
            y = which(sequences[j,]=="a" | sequences[j,]=="c" | sequences[j,]=="g" | sequences[j,]=="t")
            index = 1:length(sequences[i,])
            index = index[table(c(x,y))==2]
            seqlength=length(index)
            if(seqlength>0) dist[i,j]<-(seqlength-sum(sequences[i,index]==sequences[j,index]))/seqlength
            else{return("all sites contain missing characters")}
        }
    }
    if(nst==1) dist<--0.75*log(1-4*dist/3)
    
    return(t(dist)+dist)
}