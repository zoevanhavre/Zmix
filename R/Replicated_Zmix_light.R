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
 ReplicatedZmix_light<-function (NumRep,  sim , n, K=10,mylabels="trial", nIter=20000) {
	
 	# REPLICATE y samples
	yrep<-lapply(rep(n, NumRep),  function(x) simMe( sim, x))
 	zmixRun<-lapply(yrep,FUN=function(x, it) Zmix_lightLYRA(x, it),it=nIter)

 	docall<-do.call(rbind, lapply(zmixRun, melt))
 	K0s<-data.frame(  "Replicate"=rep(1:NumRep, each= dim(docall)[1]/NumRep)  , docall)	
	names(K0s)[-1]<-c("PT_Chain", "K0", "Proportion")
	TargetK0<-subset(K0s, PT_Chain==max(K0s$PT_Chain))

	Ymatrix<-matrix(unlist(yrep), nrow=2*NumRep, byrow=TRUE)[seq(1,2*NumRep, by=2),]  # each row is a y
	Zmatrix<-matrix(unlist(yrep), nrow=2*NumRep, byrow=TRUE)[seq(2,2*NumRep, by=2),]
# plots:
		p <- ggplot(TargetK0, aes(factor(K0), Proportion )) +  geom_boxplot(fill=NA)+xlab("Number of Groups")+ylab("Proportion of iterations")+geom_jitter(position=position_jitter(width=0.01,height=.01), alpha=.3, size=.5)+ coord_flip()
		ggsave(plot=p, filename= paste("Target_K0",mylabels,".tiff", sep="") ,
		 width=10, height=10, units='cm' )

	return(list("TargetK0"=TargetK0, "Y"=Ymatrix, "Z"=Zmatrix)) }