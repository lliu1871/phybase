run.seqgen <- function(path_seqgen="./Seq-Gen-1.3.4/seq-gen",nsim=1, seed=123, basefreq=rep(0.25,4), rate=rep(1,6), seqlength=1000, gamma=0, inv=0, treefile, saveformat="phylip", outputfile)
{
    command <- paste(path_seqgen," -mGTR ", " -z ", seed," -f", basefreq[1], ",", basefreq[2], ",", basefreq[3], ",", basefreq[4], " -r", rate[1], ",", rate[2], ",", rate[3], ",", rate[4], ",", rate[5], ",", rate[6], " -l", paste(seqlength), sep="")
    if(gamma>0) command <- paste(command, " -a", gamma, " -g 4 ", sep="")
    if(inv>0) command <- paste(command, " -i", inv, sep="")
    if(saveformat=="nexus") command <- paste(command," -on ")
    if(saveformat=="phylip") command <- paste(command," -or ")
    if(saveformat=="fasta") command <- paste(command," -of ")
    command <- paste(command, " -n", nsim," < ",  treefile, " > ", outputfile,sep="")
    print(command)
    try(system(command, intern = TRUE))
}

