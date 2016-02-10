#' Processing output of label switching function
#'
#' explain here
#' @param output of UNIVARIATE label switching function
#' @keywords postprocessing
#' @export
#' @examples
#' #nope


Zagg<-	function(USout, .Y=Y){
		
			#Pars
				.par<-melt(USout$Pars, id.vars=c("Iteration", "k"))
				theta<-aggregate( value~variable+factor(k), mean ,data=.par)
	         	K<-max(.par$k)
			
						#Pred		
				maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
	       		Zhat<- factor( apply(t(USout$Zs), 2,maxZ))
	       	 	#cbind(Zhat,USout$YZ$Y, USout$YZ$Z )   
	         	#RAND<-sum(USout$YZ$Z==Zhat)/length(Zhat)*100    
	         		# MSE
         			Zemu<-as.numeric(Zhat)
         			.Mus<-theta$value[theta$variable=="Mu"]
	         		for ( i in 1:length(Zemu)){
	         			Zemu[i]<-.Mus[as.numeric(Zhat[i])]}
	         	MSE<-sum((.Y-Zemu)^2)
	         	MAE<-sum(abs(.Y-Zemu))

	       		list("theta"=theta,  "Zpred"=Zhat, "MSE"=MSE, "MAE"=MAE)
					}
