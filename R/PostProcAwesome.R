#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope


ProstProcAWESOME<-function( Grun,  Y, prep=10000,LineUp=1,Propmin=0.3, simlabel="sim"){
		n<-length(Y$Y)  
		K<-dim(Grun$Ps)[2]	
	
		## 1. split by K0
		K0<-as.numeric(names(table(Grun$SteadyScore)))
		
		# SAVE table of tests, parameter estimates and clustering (Z's)
		p_vals<-data.frame("K0"=K0, "PropIters"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1],
			"RAND"=NA, "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
						
		GrunK0<-Grun
		#for each K0:
		for ( .K0 in 1:length(K0)){

		# split data by K0
		.iterK0<-c(1:dim(Grun$Ps)[1])[Grun$SteadyScore==K0[.K0]]
		GrunK0$Mu<-	Grun$Mu[.iterK0,]
		GrunK0$Sig<-	Grun$Sig[.iterK0,]
		GrunK0$Ps<-	Grun$Ps[.iterK0,]
		GrunK0$Loglike<-	Grun$Loglike[.iterK0]
		GrunK0$Zs<-	Grun$Zs[,.iterK0]
		GrunK0$SteadyScore<-	Grun$SteadyScore[.iterK0]

		## 2. unswitch
		GrunK0us<-QuickSwitch_allPars(GrunK0, LineUpBy=LineUp,PropMin=Propmin )

		## 3. RAND, MSE
		p_vals$RAND[.K0]<- GrunK0us$RAND
		Zetc<-Zagg(GrunK0us)
		p_vals$MAE[.K0]<- Zetc$MAE
		p_vals$MSE[.K0]<- Zetc$MSE


		## 4. Predict replicates
		PostPredFunk<-function(.GrunK0us=GrunK0us, .Zetc=Zetc){
				n<-length(Y$Y)
				swWeights<- reshape(.GrunK0us$Pars, v.names="P", idvar="Iteration", timevar="k", direction='wide', drop=c("Mu", "Sig"))[,-1]
				K<-dim(swWeights)[2]
				K<- max(as.numeric(names(table(.GrunK0us$Pars$k))))

				swMeans<- reshape(.GrunK0us$Pars, v.names="Mu", idvar="Iteration", timevar="k", direction='wide', drop=c("P", "Sig"))[,-1]
				swVariances<- reshape(.GrunK0us$Pars, v.names="Sig", idvar="Iteration", timevar="k", direction='wide', drop=c("Mu", "P"))[,-1]
			
				DrawIters<-function(x) sample(c(1:K), size=x, replace = T, prob = NULL)
				.iters<-sapply(rep(n, prep), DrawIters)
								
				# apply to .iters :   draw Z and do rnorm
				DrawRepY<-function(x){ .z<-	sample(c(1:K) ,size=1, prob=swWeights[x,]) ;cbind(rnorm(1, swMeans[x, .z], sqrt(swVariances[x, .z] )),.z )}	
				.yzrep<-sapply(.iters, DrawRepY)
				.yrep<-matrix(.yzrep[1,],nrow=prep, byrow=T)
				.zrep<-matrix(.yzrep[2,],nrow=prep, byrow=T)

				## calculate various values
				# min/max
				MinP<-sum(apply(.yrep, 1, min) < min(Y$Y))/prep
				MaxP<-sum(apply(.yrep, 1, max) > max(Y$Y))/prep

				# Prediction Concordance 
				ComputePredConcordance<-function(x){sum( (x< quantile(Y$Y, .025)) | (x > quantile(Y$Y, 1-.025))  ) /n}
				.pc<-apply(.yrep, 1, ComputePredConcordance)
				#p_vals$Concordance[.K0]<-paste(mean(.pc), " (",quantile(.pc, .025), ",", quantile(.pc, 1-.025), ")", sep="")
				Concordance<-mean(.pc)
		 		
		 		# 4.2 MSPE
				# take Z matrix and replace with estimated mean
				Zemu<-.zrep
		     	.PosteriorMeans<-.Zetc$theta$value[.Zetc$theta$variable=="Mu"]
				.PosteriorWeight<-.Zetc$theta$value[.Zetc$theta$variable=="P"]
				.PosteriorVar<-.Zetc$theta$value[.Zetc$theta$variable=="Sig"]

		         	Zemu<-apply( Zemu, c(1,2), function(x) {return(.PosteriorMeans[x])} )
			        
					MSPE_dist<-apply((.yrep-Zemu)^2, 1, sum)
					MAPE_dist<-apply(abs(.yrep-Zemu), 1, sum)

				MSPE<-mean(MSPE_dist)
				MAPE<-mean(MAPE_dist)

				### 4.3 Plot data VS replicates	
				predplot<-ggplot(data.frame("Y"=Y$Y, "Z"=Y$Z), aes(x=Y)) +theme_bw()+
				geom_line(data=melt(.yrep[,1:n]),stat="density", aes(x=value,group=Var1), size=0.5, color="blue", alpha=0.1) + geom_density(color="red", size=2)
				predplot<-predplot+ggtitle(paste("K0=", K0[.K0]))
								
				ggsave(plot=predplot, filename= paste("PredictiveDensities_",simlabel,"_K0",K0[.K0],".pdf", sep="") )
				
				return(list( "MinP"=MinP, "MaxP"=MaxP, "MAPE"=MAPE,  "MSPE"=MSPE, "Concordance"=Concordance))}

		postPredTests<-PostPredFunk( GrunK0us)
		# store output in p_vasl
		p_vals$Pmin[.K0]<-postPredTests$MinP
		p_vals$Pmax[.K0]<-postPredTests$MaxP
		p_vals$MAPE[.K0]<-postPredTests$MAPE
		p_vals$MSPE[.K0]<-postPredTests$MSPE
		p_vals$Concordance[.K0]<-postPredTests$Concordance

		}
		return(p_vals)
		}