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


getBF<-function(mydata, alpha=10, Krange=c(1:10)){
	require(bayesm)

	Data1<-list(y=matrix(mydata))
	#pick likely grid point

	km<-kmeans(Data1$y, max(Krange))
	grd<-matrix(km$centers[which.max(km$size),1])
	#grd<-matrix( mean(Data1$y))
	
	CompareLogMar<-data.frame("ncomp"=Krange,"MargDens_bm"=NA, "BF_bm"=NA  )

	for (i in 1:length(Krange)){
		ncomp<-Krange[i]
		out<-rnmixGibbs(Data=Data1,Prior=list(ncomp=ncomp, a=c(rep(alpha, ncomp))),Mcmc=list(R=20000,keep=1))
		CompareLogMar$MargDens_bm[i]<-eMixMargDen(grid=grd,probdraw= out$nmix$probdraw,compdraw= out$nmix$compdraw )	}

	for (i in 2:dim(CompareLogMar)[1])	{CompareLogMar$BF_bm[i]<- CompareLogMar$MargDens_bm[i]/CompareLogMar$MargDens_bm[i-1]}
	
	CompareLogMar2<-CompareLogMar
	for (i in 2:dim(CompareLogMar)[1])	{	
		if(CompareLogMar2$BF_bm[i]==max(na.omit(CompareLogMar2$BF_bm))){
				CompareLogMar2$BF_bm[i]<-paste(CompareLogMar2$BF_bm[i], "***")	 
		}else if(CompareLogMar2$BF_bm[i]>1){ 
				if(CompareLogMar2$BF_bm[i]!=is.na(CompareLogMar2$BF_bm[i])){
				CompareLogMar2$BF_bm[i]<-paste(CompareLogMar2$BF_bm[i], "*")	 
	} }}

	print ("                                                                     ")
	print ("                                                                     ")
	print (CompareLogMar2)
	print ("                                                                     ")
	print(" ***: maximum BF")
	print(" *: BF>1")
	print ("                                                                     ")
	niceBF<-data.frame("K"=Krange,"BF_vs_K-1"=CompareLogMar$BF_bm)[-1,]
	return(niceBF)
	}


