#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope


PostProcUnivariate_later2<-function( Grun,  mydata,LineUp=1,prep=10000,Propmin=0.05, isSim=TRUE, simlabel="sim", savelabel="PPplot", nEnd=2000){
	require(wq)
		Grun<-trimit(Out=Grun, nEnd)
		ifelse(isSim==TRUE, Y<-mydata$Y,  Y<-mydata)

		n<-length(Y)  
		K<-dim(Grun$Ps)[2]	
	
		## 1. split by K0
		K0<-as.numeric(names(table(Grun$SteadyScore)))
		
		# SAVE table of tests, parameter estimates and clustering (Z's)
		 p_vals<-data.frame("K0"=K0, "PropIters"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1],
			 "RAND"=NA, "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
						
		K0estimates<-vector("list", length(K0))
		GrunK0us_FIN<-vector("list", length(K0))

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

# PLOTS density pars
	GrunK0us$Pars$k<-as.factor(GrunK0us$Pars$k)
	p1<-ggplot(data=GrunK0us$Pars, aes(y=P, x=k)) + geom_boxplot( aes(fill=k))+ggtitle( bquote( atop(italic( .(simlabel) ), atop("Weights"))))+ ylab("")+xlab("")  +theme_bw()+  theme(legend.position = "none")
	p2<-ggplot(data=GrunK0us$Pars, aes(y=Mu, x=k)) + geom_boxplot( aes(fill=k))+ggtitle(ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means")))))+ylab("")+xlab("") +theme_bw()+  theme(legend.position = "none")
	p3<-ggplot(data=GrunK0us$Pars, aes(y=Sig, x=k)) +geom_boxplot( aes(fill=k))+ggtitle(ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), sep="")), atop("Variances")))))+ylab("")+xlab("") +theme_bw()+  theme(legend.position = "none")
	#grobframe <- arrangeGrob(p1, p2, p3, ncol=3, nrow=1,main = textGrob(paste(simlabel,": posterior parameter estimates for", K0[.K0]," groups"), gp = gpar(fontsize=8, fontface="bold.italic", fontsize=14)))
	#ggsave(plot=grobframe, filename= paste("PosteriorParDensities_",simlabel,"_K0", K0[.K0],".pdf", sep="") , width=20, height=7, units='cm' )

		ggAllocationPlot<-function( outZ, myY){
			grr<-outZ[order(myY),]
			grrTable<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)
			maxK<-max(grr)
			for (i in 1:length(myY)){rr<-factor(grr[i,], levels=1:maxK)
			grrTable<-rbind(grrTable,cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))    }
			names(grrTable)<-c("myY", "k", "Prob")
				grrTable$k<-as.factor(grrTable$k)

			gp<-ggplot(grrTable, aes(x=myY, y=k, fill=Prob)) + geom_tile()+ggtitle(  "Posterior allocations")+ 
			xlab("index of ordered y")+
			scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494" ))+theme_bw()+theme(legend.position='right')
			#ggsave( plot=gp,  filename=paste( "Allocations_", plotfilename ,"K_",maxK, ".pdf",sep="") )
			gp
			}

		p4<-ggAllocationPlot(GrunK0us$Zs, Y )
			maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
		    Zhat<- factor( apply(t(GrunK0us$Zs), 2,maxZ))
		
			## 3. RAND, MSE	
			if(isSim==TRUE){
			#maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
			#Zhat<- factor( apply(t(GrunK0us$Zs), 2,maxZ));	
			p_vals$RAND[.K0]<-(sum(mydata$Z==Zhat)/n) 
			} else { p_vals$RAND[.K0]<-'NA'}
		GrunK0us$Pars$k<-as.numeric(as.character(GrunK0us$Pars$k))

		Zetc<-ZmixUnderConstruction::Zagg(GrunK0us, Y)
		p_vals$MAE[.K0]<- Zetc$MAE
		p_vals$MSE[.K0]<- Zetc$MSE
		postPredTests<-PostPredFunk( GrunK0us,Zetc, Y, prep, simlabel)
		# store output in p_vasl
		p_vals$Pmin[.K0]<-postPredTests$MinP
		p_vals$Pmax[.K0]<-postPredTests$MaxP
		p_vals$MAPE[.K0]<-postPredTests$MAPE
		p_vals$MSPE[.K0]<-postPredTests$MSPE
		p_vals$Concordance[.K0]<-1-postPredTests$Concordance	

		p5<-postPredTests$ggp

		# CI
		.par<-melt(GrunK0us$Pars, id.vars=c("Iteration", "k"))
		theta<-aggregate( value~variable+factor(k), mean ,data=.par)
		mu<-round(aggregate( value~variable+factor(k), mean ,data=.par)[,3], 2)
		ci<-round(aggregate( value~variable+factor(k), quantile,c(0.025, 0.975) ,data=.par)[,3],2)
		thetaCI<-cbind( theta[,c(1,2)] , "value"=paste( mu, "(", ci[,1] , "," ,ci[,2] ,")", sep="" ))
		K0estimates[[.K0]]<-cbind(thetaCI, "K0"=K0[.K0])

		pdf( file= paste("PPplots_", savelabel ,"K_", K0[.K0] ,".pdf", sep=""), width=10, height=5)
 		print( wq::layOut(	list(p1, 	1, 1:2),  
	        	list(p2, 	1, 3:4),   
	         	list(p3,	1,5:6),
	         	list(p4, 	2,1:3),  
	          	list(p5, 	2,4:6)))
		dev.off()

		}
		Final_Pars<-do.call(rbind, K0estimates)
		return(list( Final_Pars, p_vals, "Z"=Zhat))
		}