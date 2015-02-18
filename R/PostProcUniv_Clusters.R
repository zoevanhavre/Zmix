#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope


PostProcUniv_Clusters<-function( Grun,  mydata,LineUp=1,Propmin=0.05, isSim=TRUE, nEnd=2000){
	require(wq)
		Grun<-trimit(Out=Grun, nEnd)
		ifelse(isSim==TRUE, Y<-mydata$Y,  Y<-mydata)

		n<-length(Y)  
		K<-dim(Grun$Ps)[2]	
	
		## 1. split by K0
		K0<-as.numeric(names(table(Grun$SteadyScore)))
		
		# SAVE table of tests, parameter estimates and clustering (Z's)
		 Props<-data.frame("K0"=K0, "PropIters"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1])
						
		K0estimates<-vector("list", length(K0))
		GrunK0us_FIN<-vector("list", length(K0))
		FinZHats<-vector("list", length(K0))
	#for each K0:
		for ( .K0 in 1:length(K0)){
		GrunK0<-Grun
		# split data by K0
		.iterK0<-c(1:dim(Grun$Ps)[1])[Grun$SteadyScore==K0[.K0]]
		GrunK0$Mu<-	Grun$Mu[.iterK0,]
		GrunK0$Sig<-	Grun$Sig[.iterK0,]
		GrunK0$Ps<-	Grun$Ps[.iterK0,]
		GrunK0$Loglike<-	Grun$Loglike[.iterK0]
		GrunK0$Zs<-	Grun$Zs[,.iterK0]
		GrunK0$SteadyScore<-	Grun$SteadyScore[.iterK0]

		## 2. unswitch
		GrunK0us<-ZmixUnderConstruction::Zswitch_univ(GrunK0, LineUp,Propmin )
		GrunK0us_FIN[[.K0]]<-GrunK0us

		maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
		 FinZHats[[.K0]]<- factor( apply(t(GrunK0us$Zs), 2,maxZ))
		
		
		}

		do.call(rbind, FinZHats)
		return(list( "ModelProbs"=Props, "Z"=FinZHats, "y"=mydata))
		}