
'test.equalgenetree'<-
function(path_raxml, inputfolder, nbootstrap)
{
    #######################################################
    # concatenation tree
    #######################################################
    allseqfile = "allsequence"
    files = list.files(inputfolder)
    reftreefiles = paste(files,".reftree",sep="")
    seqfiles = paste(inputfolder,"/",files,sep="")
    
    file.concatData(seqfiles,allseqfile)
    command = paste(path_raxml," -s",allseqfile," -mGTRGAMMA -n",allseqfile," -p",floor(runif(1)*38276153+3827262),sep="")
    try(system(command))
    
    all = read.dna.seq(allseqfile,format="phylip")
    allseq = all$seq
    allname = all$name
    alllength = dim(allseq)[2]
    alltree = read.tree(paste("RAxML_bestTree.",allseqfile,sep=""))
    ngene = length(seqfiles)
    seqlength = 1:ngene
    
    #######################################################
    # fit data to the concatenation tree
    #######################################################
    for(i in 1:ngene)
    {
        y<-read.dna.seq(seqfiles[i],format="phylip")
        seqlength[i] = dim(y$seq)[2]
        genename<-y$name
        delname<-allname[which(is.na(match(allname,genename))==1)]
        newtree<-write.tree(drop.tip(alltree,delname))
        write(newtree,reftreefiles[i])
        
        command = paste(path_raxml," -s",seqfiles[i]," -mGTRGAMMA -nfix",files[i], " -f e -t",reftreefiles[i]," -p",floor(runif(1)*38276153+3827262),sep="")
        try(system(command))
        command = paste(path_raxml," -s", seqfiles[i]," -mGTRGAMMA -n",files[i]," -p",floor(runif(1)*38276153+3827262),sep="")
        try(system(command))
    }
    
    #######################################################
    # test statistic
    #######################################################
    fixloglike=loglike=1:ngene
    for(i in 1:ngene)
    {
        x<-scan(paste("RAxML_log.",files[i],sep=""),what="character")
        loglike[i]<-as.numeric(x[length(x)])
        x<-scan(paste("RAxML_log.fix",files[i],sep=""),what="character")
        fixloglike[i]<-as.numeric(x[length(x)])
    }
    teststat=sum(loglike-fixloglike)
    teststat_gene = loglike-fixloglike
    
    #######################################################
    # null distribution of teststat using bootstrap
    #######################################################
    currentfolder=getwd()
    for(i in 1:nbootstrap)
    {
        try(system(paste("mkdir",i)))
        for(j in 1:ngene)
        {
            y<-sort(read.tree.string(reftreefiles[j],format="phylip")$name)
            index<-match(y,allname)
            s<-sample(1:alllength,size=seqlength[j],replace=TRUE)
            write.dna.seq(allseq[index,s], allname[index], paste("./",i,"/",files[j],sep=""),format="phylip")
        }
        setwd(paste("./",i,sep=""))
	x=list.files()
        file.concatData(x, allseqfile)
        command = paste(path_raxml," -s",allseqfile," -mGTRGAMMA -n",allseqfile," -p",floor(runif(1)*38276153+3827262),sep="")
        try(system(command))
        alltree1 = read.tree(paste("RAxML_bestTree.",allseqfile,sep=""))
        
        for(j in 1:ngene)
        {
            y<-read.dna.seq(files[j],format="phylip")
            genename<-y$name
            delname<-allname[which(is.na(match(allname,genename))==1)]
            newtree<-write.tree(drop.tip(alltree1,delname))
            write(newtree,reftreefiles[j])
            
            command = paste(path_raxml," -s",files[j]," -mGTRGAMMA -nfix",files[j], " -f e -t",reftreefiles[j]," -p",floor(runif(1)*38276153+3827262),sep="")
            try(system(command))
            command = paste(path_raxml," -s", files[j]," -mGTRGAMMA -n",files[j]," -p",floor(runif(1)*38276153+3827262),sep="")
            try(system(command))
        }
        setwd(currentfolder)
        
    }
    
    #######################################################
    # summarize result
    #######################################################
    folder<-paste("./",1:nbootstrap,sep="")
    currentfolder<-getwd()
    file<-paste("RAxML_log.",files,sep="")
    fixfile<-paste("RAxML_log.fix",files,sep="")
    bootloglike<-1:ngene
    bootfixloglike<-1:ngene
    
    test_bootstrap<-matrix(0,nbootstrap,ngene)
    for(j in 1:nbootstrap)
    {
        setwd(folder[j])
        for(i in 1:ngene)
        {
            x<-scan(file[i],what="numeric")
            bootloglike[i]<-as.numeric(x[length(x)])
            
            x<-scan(fixfile[i],what="numeric")
            bootfixloglike[i]<-as.numeric(x[length(x)])
        }
        test_bootstrap[j,]<-bootloglike-bootfixloglike
        setwd(currentfolder)
    }
    
    testb = apply(test_bootstrap,1,sum)
    pvalue = sum(testb > teststat)/nbootstrap
    pvalue.gene = 1:ngene
    for(i in 1:length(teststat_gene))
	pvalue.gene[i] = sum(test_bootstrap[,i] > teststat_gene[i])/nbootstrap
    
    z = list(teststat=teststat, testb=testb, pvalue=pvalue, pvalue.gene=pvalue.gene)
    z
}

