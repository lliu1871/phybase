'run.modeltest' <-
function(path_jmodeltest="./jmodeltest2/dist/jModelTest.jar", seqfile, nmodel=3, outputfile)
{
    command <- paste("java -jar ", path_jmodeltest," -d", seqfile, " -s", nmodel, " -g 4 -i -f -AIC >", outputfile)
    result<-try(system(command))
}

