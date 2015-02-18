#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope


Zmix_switchByK0<-function( Grun,  mydata, prep=10000,LineUp=1,Propmin=0.05, isSim=TRUE, simlabel="sim"){
		
		ifelse(isSim==TRUE, Y<-mydata$Y,  Y<-mydata)

		n<-length(Y)  
		K<-dim(Grun$Ps)[2]	
	
		## 1. split by K0
		K0<-as.numeric(names(table(Grun$SteadyScore)))
		
		# SAVE table of tests, parameter estimates and clustering (Z's)
		p_vals<-data.frame("K0"=K0, "PropIters"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1],
			"RAND"=NA, "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
						
		K0estimates<-vector("list", length(K0))
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
		GrunK0us<-Zswitch_univ(GrunK0, LineUp,Propmin )

# PLOTS
p1<-ggplot(data=GrunK0us$Pars, aes(x=P, fill=factor(k))) + geom_density( alpha=0.4)+ggtitle("Weights ")+ylab("")+xlab("")  +  theme(legend.position = "none")
p2<-ggplot(data=GrunK0us$Pars, aes(x=Mu, fill=factor(k))) + geom_density( alpha=0.4)+ggtitle("Means")+ylab("")+xlab("") +  theme(legend.position = "none")
p3<-ggplot(data=GrunK0us$Pars, aes(x=Sig, fill=factor(k))) +geom_density(alpha=0.4)+ggtitle("Variance")+ylab("")+xlab("") +  theme(legend.position = "none")
grobframe <- arrangeGrob(p1, p2, p3, ncol=3, nrow=1,main = textGrob(paste(simlabel,": posterior parameter estimates for", K0[.K0]," groups"), gp = gpar(fontsize=8, fontface="bold.italic", fontsize=14)))
ggsave(plot=grobframe, filename= paste("PosteriorParDensities_",simlabel,"_K0", K0[.K0],".tiff", sep="") , width=20, height=7, units='cm' , compression="lzw")

		## 3. RAND, MSE	
		if(isSim==TRUE){
			maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
			Zhat<- factor( apply(t(GrunK0us$Zs), 2,maxZ))
			p_vals$RAND[.K0]<-(sum(mydata$Z==Zhat)/n)*100
						} else { p_vals$RAND[.K0]<-'NA'}

		Zetc<-ZmixUnderConstruction::Zagg(GrunK0us, Y)
		p_vals$MAE[.K0]<- Zetc$MAE
		p_vals$MSE[.K0]<- Zetc$MSE

		K0estimates[[.K0]]<-cbind(Zetc$theta, "K0"=K0[.K0])

		## 4. Predict replicates
		
		postPredTests<-ZmixUnderConstruction::PostPredFunk( GrunK0us,Zetc, Y, prep, simlabel)
		# store output in p_vasl
		p_vals$Pmin[.K0]<-postPredTests$MinP
		p_vals$Pmax[.K0]<-postPredTests$MaxP
		p_vals$MAPE[.K0]<-postPredTests$MAPE
		p_vals$MSPE[.K0]<-postPredTests$MSPE
		p_vals$Concordance[.K0]<-1-postPredTests$Concordance		}
		return(list(p_vals, K0estimates, GrunK0us))
		}