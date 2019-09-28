'parse.raxml' <-
function(raxml_info_file="RAxML_info.out")
{
    x<-scan(raxml_info_file,what="character",sep="\n")
    y<-x[grep("Final GAMMA-based Score of best tree",x)]
    loglike <- as.numeric(gsub("Final GAMMA-based Score of best tree","",y))
    basefreq<-gsub("Base frequencies:","",x[grep("Base frequencies",x)])
    b=as.numeric(unlist(strsplit(basefreq,split=" ")))
    basefreq = b[is.na(b)==FALSE]
    y<-x[grep("alpha",x)]
    y<-unlist(strsplit(y[2],split=" "))
    alpha<-y[2]
    rates<-y[10:15]
    
    #multiple runs
    bestrun=as.numeric(unlist(strsplit(gsub("Starting final GAMMA-based thorough Optimization on tree ","",x[grep("Starting final GAMMA-based thorough Optimization on tree",x)]),split=" "))[1])
    if(length(bestrun)>0){
        y=x[grep(paste("Inference\\[",bestrun,"\\]",sep=""),x)[1]+1]
        y=gsub("rates\\[0\\] ac ag at cg ct gt: ","",gsub("alpha\\[0\\]: ","",y))
        y=as.numeric(unlist(strsplit(y,split=" ")))
        alpha = y[1]
        rates = y[2:7]
    }
    tree<-read.tree.string(gsub("info","bestTree",raxml_info_file),format="phylip")
    
    z=list(basefreq=basefreq, alpha=alpha, rates=rates, tree=tree$tree, loglike=loglike)
    z
}

