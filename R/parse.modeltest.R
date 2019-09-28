'parse.modeltest' <-
function(outputfile)
{
    result<-scan(outputfile,what="character",sep="\n")
    x <- result[grep("Best Models", result)+4]
    x <- unlist(strsplit(x,split="\t"))
    x <- x[x!=""]
    
    model <- x[2]
    basefreq <- as.numeric(x[3:6])
    kappa <- as.numeric(x[7])
    tstv <- as.numeric(x[8])
    ratematrix <- as.numeric(x[9:14])
    gamma <- 0
    inp <- 0
    if(x[15] != "N/A") inp <- as.numeric(x[15])
    if(x[16] != "N/A") gamma <- as.numeric(x[16])
    
    loglike <- -as.numeric(unlist(strsplit(result[grep("Model selected",result)+3],split="="))[2])
    tree <- unlist(strsplit(result[grep("for the best AIC model ",result)],split="="))[2]
    
    x <- grep("Selection uncertainty",result)+4
    y <- grep("negative log likelihod",result)-2
    z <- unlist(strsplit(result[x-2],split=" "))
    name = z[z!=""]
    outtable <- matrix("",nrow = y-x+1, ncol=length(name))
    colnames(outtable) <- name
    
    for(i in x:y)
    {
        z<-unlist(strsplit(result[i],split=" "))
        outtable[i-x+1,] <- z[z!=""]
    }
    
    m <- list(model=as.character, loglike=as.numeric, basefreq=as.numeric, kappa=as.numeric, tstv=as.numeric, ratematrix=as.numeric, gamma=as.numeric, inp=as.numeric, tree=as.character, table=as.matrix)
    m$model <- model
    m$loglike <- loglike
    m$basefreq <- basefreq
    m$kappa <- kappa
    m$tstv <- tstv
    m$ratematrix <- ratematrix
    m$gamma <- gamma
    m$inp <- inp	
    m$tree <- tree
    m$table <- outtable
    
    m
}