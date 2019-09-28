'test.submodel.valid.raxml' <-
function(inputfile, path_raxml, path_seqgen, nbootstrap=100)
{
    simfile<-paste(inputfile,".sim",sep="")
    command <- paste(path_raxml,"-N 10 -m GTRGAMMA -n",inputfile, " -s", inputfile," -mGTRGAMMA","-p", floor(runif(1)*5945478+3843623))
    try(system(command))
    fit = parse.raxml(paste("RAxML_info.",inputfile,sep=""))
    
    run.seqgen(path_seqgen,
    seed = floor(4542354*runif(1)+45439587), basefreq = fit$basefreq, rate = fit$rates, seqlength = "100000",
    gamma = fit$alpha, treefile=paste("RAxML_bestTree.",inputfile,sep=""), saveformat = "phylip", outputfile=simfile)
    
    
    simdata<-read.dna.seq(simfile,format="phylip")
    name<-sort(simdata$name)
    nspecies<-length(name)
    simseq<-simdata$seq
    simseq<-tolower(simseq)
    simseq<-simseq[match(name,simdata$name),]
    simlength<-dim(simseq)[2]
    
    data<-read.dna.seq(inputfile,format="phylip")
    seq<-tolower(data$seq)
    seq<-seq[match(name,data$name),]
    seqlength<-dim(seq)[2]
    
    ######################################################
    # chisquare-test for base frequencies
    ######################################################
    base.test=matrix(0,nrow=nspecies,ncol=11)
    colnames(base.test) = c("species","obs_a","obs_c","obs_g","obs_t","exp_a","exp_c","exp_g","exp_t","teststat","p.value")
    base.test[,1] = name
    for(i in 1:nspecies)
    {
        obs<-c(sum(seq[i,]=="a"),sum(seq[i,]=="c"),sum(seq[i,]=="g"),sum(seq[i,]=="t"))
        base.test[i,2:5] <- obs
        base.test[i,6:9] <- fit$basefreq
        result<-chisq.test(obs,p=fit$basefreq,rescale.p =T)
        base.test[i,10] <- round(result$statistic,5)
        base.test[i,11] <- round(result$p.value,5)
    }
    
    print("The Marginal Test for base frequencies is complete!")
    ######################################################
    # chisquare-test for doublet frequencies
    ######################################################
    pair_pvalue<-matrix(0,nspecies,nspecies)
    colnames(pair_pvalue) = name
    rownames(pair_pvalue) = name
    for(m in 1:(nspecies-1))
    for(n in (m+1):nspecies)
    {
        simseq1<-simseq[m,]
        simseq2<-simseq[n,]
        simpair<-.pairfreq(simseq1,simseq2)
        exppair<-simpair/sum(simpair)
        obsseq1<-seq[m,]
        obsseq2<-seq[n,]
        obspair<-.pairfreq(obsseq1,obsseq2)
        if(length(which(exppair==0))>0) obspair = obspair[-which(exppair==0)]
        if(length(which(exppair==0))>0) exppair = exppair[-which(exppair==0)]
        if(sum(obspair)==0){pair_pvalue[m,n]=NA;next;}
        pair_pvalue[m,n]<-chisq.test(obspair,p=exppair,rescale.p=TRUE, simulate.p.value=T)$p.value
    }
    
    print("The Marginal Test for doublet frequencies is complete!")
    ######################################################
    # chisquare-test for site patterns
    ######################################################
    newseq <- gsub("N", "?", seq)
    newseq <- gsub("-", "?", newseq)
    newseq <- gsub("O", "?", newseq)
    site_pattern <- apply(newseq, 2, function(x) paste(x, collapse = "", sep = ""))
    obs_freq <- table(site_pattern)
    pattern <- names(obs_freq)
    
    new = .patternfreq(pattern, newseq)
    obs_pattern_freq <- new$freq
    obs_pattern <- new$pattern
    exp_pattern_freq = .patternfreq(obs_pattern, simseq)$freq
    teststat<-sum(abs(exp_pattern_freq-obs_pattern_freq))
    
    nbootstrap<-nbootstrap
    bootstat<-1:nbootstrap
    for(j in 1:nbootstrap)
    {
        bootsample<-simseq[,sample(1:simlength,seqlength)]
        freq = .patternfreq(obs_pattern, bootsample)$freq
        bootstat[j]<-sum(abs(exp_pattern_freq-freq))
    }
    sitepattern.pvalue = sum(bootstat>teststat)/nbootstrap
    
    z <- list(base.test=base.test, doublet.pvalue=pair_pvalue, sitepattern.pvalue=sitepattern.pvalue)
    z
}


.patternfreq<-function(pattern, seq)
{
    nspecies = dim(seq)[1]
    newobs_freq = 1:length(pattern)
    names(newobs_freq) = pattern
    for(i in 1:length(pattern))
    {
        str<-unlist(strsplit(pattern[i],split=""))
        index<-which(str=="?")
	if(length(index)==nspecies){
	    newobs_freq[i] <- NA
	    next
	}
        {if(length(index)>0){
            newp<-paste(str[-index],collapse="",sep="")
            news<-seq[-index,]
        }else{
            newp<-paste(str,collapse="",sep="")
            news<-seq
        }}
        if(length(index)<nspecies-1) news<-apply(news, 2, function(x) paste(x, collapse = "", sep = ""))
        nmissing<-length(grep("\\?",news))
        {if(nmissing<length(news)*0.8) newobs_freq[i]<-length(grep(newp,news))/(length(news)-nmissing)
        else newobs_freq[i] <- NA}
    }
    newobs<-newobs_freq[!is.na(newobs_freq)]
    newpattern<-names(newobs)
    
    z <- list(pattern=newpattern, freq=newobs)
    z
}

.pairfreq<-function(seq1,seq2)
{
    x<-paste(seq1,seq2,sep="")
    result<-c(sum(x=="aa"),sum(x=="ac"),sum(x=="ag"),sum(x=="at"),sum(x=="ca"),sum(x=="cc"),sum(x=="cg"),sum(x=="ct"),sum(x=="ga"),sum(x=="gc"),sum(x=="gg"),sum(x=="gt"),sum(x=="ta"),sum(x=="tc"),sum(x=="tg"),sum(x=="tt"))
    return(result)
}

