'parse.phyml' <-
function(phyml_stats_file="phyml_stats.txt")
{
    x<-scan(phyml_stats_file,what="character",sep="\n")
    y<-x[grep("Log-likelihood",x)]
    loglike <- as.numeric(gsub(". Log-likelihood:","",y))
    basefreq<-1:4
    basefreq[1]<-as.numeric(gsub("- f\\(A\\)=","",x[grep("- f\\(A\\)=",x)]))
    basefreq[2]<-as.numeric(gsub("- f\\(C\\)=","",x[grep("- f\\(C\\)=",x)]))
    basefreq[3]<-as.numeric(gsub("- f\\(G\\)=","",x[grep("- f\\(G\\)=",x)]))
    basefreq[4]<-as.numeric(gsub("- f\\(T\\)=","",x[grep("- f\\(T\\)=",x)]))
    alpha<-as.numeric(gsub("- Gamma shape parameter:","",x[grep("- Gamma shape parameter:",x)]))
    inv<-as.numeric(gsub(". Proportion of invariant:","",x[grep("Proportion of invariant:",x)]))
    tv<-as.numeric(gsub(". Transition/transversion ratio:","",x[grep("transversion ratio:",x)]))
    if(length(tv)==0)
    {
    rates = 1:6
    rates[1]<-as.numeric(gsub("A <-> C","",x[grep("A <-> C",x)]))
    rates[2]<-as.numeric(gsub("A <-> G","",x[grep("A <-> G",x)]))
    rates[3]<-as.numeric(gsub("A <-> T","",x[grep("A <-> T",x)]))
    rates[4]<-as.numeric(gsub("C <-> G","",x[grep("C <-> G",x)]))
    rates[5]<-as.numeric(gsub("C <-> T","",x[grep("C <-> T",x)]))
    rates[6]<-as.numeric(gsub("G <-> T","",x[grep("G <-> T",x)]))
    }
    tree<-read.tree.string(gsub("stats","tree",phyml_stats_file),format="phylip")
    
    z=list(basefreq=basefreq, alpha=alpha, inv=inv, tv=tv, rates=rates, tree=tree$tree, loglike=loglike)
    z
}

