#' Bayes Factors Using Bayesm 
#'
#' This function draws samples from a Wishart dist
#' @param mydata, alpha, Krange
#' @keywords bayes factor, bayesm
#' @export
#' @examples
#' # see bayesm documentation
#' sim1_n100<-sim1func(n=300) 
#' BF<-getBF(sim1_n100$Y, 1, c(1:6))
#' BF


getBFk0<-function(mydata, runlength){
	 alpha<-1; Krange<-c(1:10)
	require(bayesm)
	Data1<-list(y=matrix(mydata$Y))
	#pick likely grid point
	km<-kmeans(Data1$y, max(Krange))
	grd<-matrix(km$centers[which.max(km$size),1])
	
	
	CompareLogMar<-data.frame("ncomp"=Krange,"MargDens_bm"=NA, "BF_bm"=NA  )

	for (i in 1:length(Krange)){
		ncomp<-Krange[i]
		out<-rnmixGibbs(Data=Data1,Prior=list(ncomp=ncomp, a=c(rep(alpha, ncomp))),Mcmc=list(R=runlength,keep=1))
		CompareLogMar$MargDens_bm[i]<-eMixMargDen(grid=grd,probdraw= out$nmix$probdraw,compdraw= out$nmix$compdraw )	}

	for (i in 2:dim(CompareLogMar)[1])	{CompareLogMar$BF_bm[i]<- CompareLogMar$MargDens_bm[i]/CompareLogMar$MargDens_bm[i-1]}
	
	niceBF<-data.frame("K"=Krange,"BF_vs_K-1"=CompareLogMar$BF_bm)[-1,]

	BFK0<-1
	if (max(niceBF[,2]>1)){ BFK0<- niceBF[which.max(niceBF[,2]) ,1]}

	return(BFK0)
	}


