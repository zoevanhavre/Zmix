 #' A replicated  Function
#'
#' ...
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))
#' 
#'   
#'  go<- ReplicatedZmix_light(3, 1, 20, K=5, iter=400)
#'
#'
#'
 ReplicatedZmix_parallelQUICK<-function (NumRep,  sim , n, mylabels="trial", ...) {
	require('parallel')
 	# REPLICATE y samples
 	dtc<-min(detectCores(), NumRep)
 	print(paste ( "Using ",dtc, " Cores"))
	yrep<-lapply(rep(n, NumRep),  function(x) simMe( sim, x))
 	#zmixRun<-lapply(yrep, function(x){   Zmix_lightLYRA(x, K,...)} )
 	zmixRun<-mclapply(yrep, FUN= function(x) Zmix_lightLYRAquicktry(x) , mc.cores=dtc) 
 	
 	docall<-do.call(rbind, lapply(zmixRun, melt))
 	K0s<-data.frame(  "Replicate"=rep(1:NumRep, each= dim(docall)[1]/NumRep)  , docall)	
	names(K0s)[-1]<-c("PT_Chain", "K0", "Proportion")
	TargetK0<-subset(K0s, PT_Chain==max(K0s$PT_Chain))
	
		# plots:
		p <- ggplot(TargetK0, aes(factor(K0), Proportion )) +  geom_boxplot()+xlab("Number of Groups")+ylab("Proportion of iterations")+geom_jitter(position=position_jitter(width=0.01,height=.01), alpha=.3, size=.5)+ coord_flip()
		ggsave(plot=p, filename= paste("Target_K0",mylabels,".tiff", sep="") ,
		 width=10, height=10, units='cm' )

	Ymatrix<-matrix(unlist(yrep), nrow=2*NumRep, byrow=TRUE)[seq(1,2*NumRep, by=2),]  # each row is a y
	Zmatrix<-matrix(unlist(yrep), nrow=2*NumRep, byrow=TRUE)[seq(2,2*NumRep, by=2),]
	save(TargetK0, file="ifIexistthisworks.RDATA")
	return(list("TargetK0"=TargetK0, "K0s"=K0s, "Y"=Ymatrix, "Z"=Zmatrix )) }
	