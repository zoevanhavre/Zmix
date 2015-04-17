#' Process the output of Zmix_univ_tempered
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope


Process_Output_Zmix<-function( Grun, isSim=FALSE,  Burn=2000, LineUp=1, Pred_Reps=100, Zswitch_Sensitivity=0.01, makePlots=TRUE, Plot_Title="Results", SaveFileName="zmix", PlotType="Boxplot"){
		mydata<-Grun$YZ
		require(wq)
		Grun<-trimit(Out=Grun, Burn)
		ifelse(isSim==TRUE, Y<-mydata$Y,  Y<-mydata)

		n<-length(Y)  
		K<-dim(Grun$Ps)[2]	
	minq4PLOT<-0.01
    	both<-cbind(melt(Grun$Mu)[,3], log(melt(Grun$Sig)[,3]^2))
    	Weights<-melt(Grun$Ps)[,3]
	color <- as.factor(melt(Grun$Mu)[,2])
	raincol<-rainbow(length(levels(color))); 	levels(color)<-raincol
	trancol<-sapply(c(1:length(color)), function(x) adjustcolor(color[x], alpha.f=Weights[x]))
	minmaxMEANS<-c(min(both[Weights>minq4PLOT ,1]), max(both[Weights>minq4PLOT ,1]))
# PLOT 1
slices <- prop.table(table((factor(Grun$SteadyScore, levels=c(1:K)))))
Lab.palette <- colorRampPalette(rainbow(K*3, alpha=.3), space = "Lab")

if(makePlots==TRUE){
pdf(file=paste(SaveFileName, "_MCMCpp.pdf",sep='') ,width=8, height=3)
	par(mfrow=c(1, 3))
	barplot(slices, ylim=c(0,1), main="Number of non-empty states", xlab="Number of non-empty states", ylab="Probability (from MCMC)") 
	abline(h=seq(0, 1, .05), lwd=0.5, col='LightGrey')
	smoothScatter(both, colramp = Lab.palette, main="Posterior Surface", nrpoints = 0, xlab="Mean", ylab="LOG(Variance)")
	plot(both, col=trancol, xlim=minmaxMEANS, ylim=c(0,100), xlab="Mean", ylab="Variance", bg='grey', main=paste("Posterior Samples:" ,"\n Transparency By Weight"))
	#if(is.na(trueValues)==FALSE){	points(trueValues, pch=7, cex=2)}
dev.off()
}


		## 1. split by number of components
		K0<-as.numeric(names(table(Grun$SteadyScore)))

		# SAVE table of tests, parameter estimates and clustering (Z's)
		 p_vals<-data.frame("K0"=K0, "Probability"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1],
			"MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
						
		K0estimates<-vector("list", length(K0))
		Zestimates<-vector("list", length(K0))
		GrunK0us_FIN<-vector("list", length(K0))
		ZTable<-vector("list", length(K0))
	#for each K0:
		for ( .K0 in 1:length(K0)){
	if( p_vals$Probability[.K0]>0.05){
		GrunK0<-Grun
		# split data by K0
		.iterK0<-c(1:dim(Grun$Ps)[1])[Grun$SteadyScore==K0[.K0]]
		GrunK0$Mu<-	Grun$Mu[.iterK0,]
		GrunK0$Sig<-Grun$Sig[.iterK0,]
		GrunK0$Ps<-	Grun$Ps[.iterK0,]
		GrunK0$Loglike<-Grun$Loglike[.iterK0]
		GrunK0$Zs<-	Grun$Zs[,.iterK0]
		GrunK0$SteadyScore<-Grun$SteadyScore[.iterK0]

		## 2. unswitch
		GrunK0us<-Zswitch(GrunK0, LineUp, Zswitch_Sensitivity )
		GrunK0us_FIN[[.K0]]<-GrunK0us

# PLOTS density pars
if(makePlots==TRUE){
	GrunK0us$Pars$k<-as.factor(GrunK0us$Pars$k)
	if(PlotType=='Density'){
	p1<-ggplot(data=GrunK0us$Pars, aes(x=P, fill=k)) + geom_density( alpha=0.4)+ggtitle( bquote( atop(italic( .(Plot_Title) ), atop("Weights"))))+ ylab("")+xlab("")  +theme_bw()+  theme(legend.position = "none")
	p2<-ggplot(data=GrunK0us$Pars, aes(x=Mu, fill=k)) + geom_density( alpha=0.4)+ggtitle(ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means")))))+ylab("")+xlab("") +theme_bw()+  theme(legend.position = "none")
	p3<-ggplot(data=GrunK0us$Pars, aes(x=Sig, fill=k)) +geom_density(alpha=0.4)+ggtitle(ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=", .(p_vals$Probability[.K0]), sep="")), atop("Variances")))))+ylab("")+xlab("") +theme_bw()+  theme(legend.position = "none")
	#grobframe <- arrangeGrob(p1, p2, p3, ncol=3, nrow=1,main = textGrob(paste(Plot_Title,": posterior parameter estimates for", K0[.K0]," groups"), gp = gpar(fontsize=8, fontface="bold.italic", fontsize=14)))
	#ggsave(plot=grobframe, filename= paste("PosteriorParDensities_",Plot_Title,"_K0", K0[.K0],".pdf", sep="") , width=20, height=7, units='cm' )
} else if (PlotType=="Boxplot"){
	pii.mean = aggregate(P ~ k, GrunK0us$Pars, mean)
	mu.mean = aggregate(Mu ~ k, GrunK0us$Pars, mean)
	var.mean = aggregate(Sig ~ k, GrunK0us$Pars, mean)

		p1<-ggplot(data=GrunK0us$Pars, aes(y=P, x=k)) + geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle( bquote( atop(italic( .(Plot_Title) ), atop("Weights"))))#+ geom_text(data =pii.mean, aes(label=signif(P,4)),size=4,  col='yellow',vjust = 1)
		p2<-ggplot(data=GrunK0us$Pars, aes(y=Mu, x=k))+ geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle(ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means")))))
		p3<-ggplot(data=GrunK0us$Pars, aes(y=Sig, x=k)) + geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle(ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=", .(p_vals$Probability[.K0]), sep="")), atop("sqrt(Variance).")))))
}
	}

	# ALOC  PROBABILITIES
Ztemp<-GrunK0us$Zs

ZTable[[.K0]]<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)
			maxK<-max(Ztemp)
			for (i in 1:dim(Ztemp)[1]){rr<-factor(Ztemp[i,], levels=1:maxK)
			ZTable[[.K0]]<-rbind(ZTable[[.K0]],cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))    }
			names(ZTable[[.K0]])<-c("Yid", "k", "Prob")
				ZTable[[.K0]]$k<-as.factor(ZTable[[.K0]]$k)


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

		if(makePlots==TRUE){ p4<-ggAllocationPlot(GrunK0us$Zs, Y )}
		
		maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
		    Zhat<- factor( apply(t(GrunK0us$Zs), 2,maxZ))
		  Zestimates[[.K0]]<-Zhat
			## 3. , MSE	
			
		GrunK0us$Pars$k<-as.numeric(as.character(GrunK0us$Pars$k))

		Zetc<-ZmixUnderConstruction::Zagg(GrunK0us, Y)
		p_vals$MAE[.K0]<- Zetc$MAE
		p_vals$MSE[.K0]<- Zetc$MSE
		postPredTests<-PostPredFunk( GrunK0us,Zetc, Y, Pred_Reps, Plot_Title)
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
if(makePlots==TRUE){
if(K0[.K0]>1){
		
		pdf( file= paste("PPplots_", SaveFileName ,"K_", K0[.K0] ,".pdf", sep=""), width=10, height=5)
 		print( wq::layOut(	list(p1, 	1, 1:2),  
	        	list(p2, 	1, 3:4),   
	         	list(p3,	1,5:6),
	         	list(p4, 	2,1:3),  
	          	list(p5, 	2,4:6)))
		dev.off()
		}}

		}}
		Final_Pars<-do.call(rbind, K0estimates)
		print(p_vals)
		#Result<-list( Final_Pars, p_vals, "Z"=Zhat)
	#save(Result, file=paste("PPresults_", SaveFileName ,".RDATA", sep=""))
		return(list( Final_Pars, p_vals, Zestimates, ZTable, "Pars_us"=GrunK0us_FIN))
		}
