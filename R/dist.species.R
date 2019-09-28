'dist.species'<-
function (dist, species.structure)
{
	nspecies<-dim(species.structure)[1]
   	dis<-matrix(0,nspecies,nspecies)
   	for (i in 1:(nspecies-1)) {
        	for(j in (i+1):nspecies){
			dis[i,j]<-mean(dist[which(species.structure[i,]==1), which(species.structure[j,]==1)], na.rm=TRUE)

   	}}
	dis<-dis+t(dis)
    	dis
}

