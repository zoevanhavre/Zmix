#' Label switch fix for Univariate Normals 
#'
#' This function fixes the label switching in the output of a univariate Gibbs sampler. Requires a constant number of non-empty groups in samples, but can deal with empty groups.
#' @param out_trim: LIST of output of Gibbs sampler with the following components:
#' 
#' 		[[Zs]] estimated component (allocation) for each Y at each iteration (rows=n by columns=Iterations)  
#'
#' 		[[Mu]] matrix of mu samples  (rows=iterations, columns  = components), 
#'
#'		[[Sig]] matrix variance  (really Sig^2, but doesnt matter here anyway)
#'
#' 		[Ps]] matrix weights
#'
#' 		[[Loglike]] Loglikelihood (vector the length of iterations). Change if needed, currently used to find reference iteration (max(loglikelihood)).  
#'
#' 		[[SteadyScore]] vector of number of non-empty groups at each iteration  (Not important, ok if not there)
#'
#'
#' 		
#' @param LineUpBy : currently can be set to 1 or 2. It is only used to rename the FINAL output 
#' 		1: final output will order Components by posterior weight (Z=1 for max(weight),etc)
#' 		2: final output ordered by variance of components
#' 	Again, note this doesnt affect the actualy unswitching, which uses no ordering. it simply relabels the reference groups according to chosen parameter to nice output.
#'
#'
#' @param PropMin = Minimum allowable proportion of a group required to activate second level check., 
#' 		Usually 0.1 sufficient.
#' 		Smaller values increase time function takes, but are needed for Smaller samples or very close components (try 0.05)
#' 	If function gives error, this is first thing to try (make it smaller)
#'
#' @keywords label switching univariate gaussian gibbs
#'
#' @export
#'
#' @examples
#' #
#' # if you have an unsteady number of posterior  non-empty groups:
#' # subset posterior samples first by this then unswitch.
#'
#' # Simple example:
#'
#'
#' set.seed(88); library(ggplot2)
#' # Make simple mixture dataset, 
#' 	dat1<-c(rnorm(50, mean=1), rnorm(100, mean=5))
#' # Run gibbs sampler of choice 
#'  run1<-Zmix_univ_tempered(dat1,isSim=FALSE,  iter=1000, k=10, alphas= c(30, 20, 10, 5, 3, 1, 0.5, 1/2^(c(2,3,4,5,6, 8, 10, 15, 20, 30))))
#'	run1<- trimit(run1, nEnd=500) # trim
#'	K0<-as.numeric(names(table(run1$SteadyScore))) # check out number of non-empty groups in iterations, if available
#' 	
#'  # Undo the LABEL SWITCHING
#' 	runUS<-QuickSwitch_allPars(run1, LineUpBy=1,PropMin=0.1 )
#'
#'  # Check out results
#'	p1<-ggplot(data=runUS$Pars, aes(x=Iteration, y=P, group=factor(k), colour=factor(k))) + geom_line() + geom_point()+ggtitle("Unswitched Weights")
#'	p2<-ggplot(data=runUS$Pars, aes(x=Iteration, y=Mu, group=factor(k), colour=factor(k))) + geom_line() + geom_point()+ggtitle("Means")
#'	p3<-ggplot(data=runUS$Pars, aes(x=Iteration, y=Sig, group=factor(k), colour=factor(k))) + geom_line() + geom_point()+ggtitle("Variance")
#' 	multiplot(p1,p2,p3)





	QuickSwitch_allPars<-function(out_trim, LineUpBy=1,PropMin=0.1 ){
			K<-dim(out_trim$Ps)[2]
			
			#ifelse(isSim==TRUE, Y<-mydata$Y, Y<-mydata)
			
			# Pick Reference = Max log Likelihood
			wml<-which.max(out_trim$Loglike)
			Zref<-factor(out_trim$Zs[,	wml], levels=1:K,ordered=FALSE) 
			if (LineUpBy==1){	
				FinalOrderChoice<-order(out_trim$Ps[wml,], decreasing=TRUE)		
			#	refPar<-out_trim$Ps[wml,FinalOrderChoice]##***##
				# comparing pars:
				non0ref<-FinalOrderChoice[1:sum(table(Zref)>0)]
				refComp<-c(out_trim$P[wml,non0ref], out_trim$Mu[wml,non0ref], out_trim$Sig[wml,non0ref])

			} else if(LineUpBy==2){	
				.tbs<-table(Zref)
				.tbs[.tbs>0]<-1
				FinalOrderChoice	 <-order(.tbs*out_trim$Sig[wml,], decreasing=TRUE)
			# PICK SUPPORTING PARS
			#refPar<-out_trim$Sig[wml,FinalOrderChoice]##***##

			# comparing pars:
			non0ref<-FinalOrderChoice[1: sum(.tbs)]  # not right
			refComp<-c(out_trim$P[wml,non0ref], out_trim$Mu[wml,non0ref], out_trim$Sig[wml,non0ref])

			}
			#levels(Zref)<-FinalOrderChoice
			levels(Zref)<- c(1:K)[order(FinalOrderChoice)]
			Zref<- factor(Zref,levels(Zref)[order(levels(Zref))])
			
			
			# storage dataframes:
			AllPars<-data.frame('Iteration'=NA, 'k'=NA, 'Ps'=NA, 'Mu'=NA, 'Sig'=NA)[numeric(0), ]
			Zfixed<-out_trim$Zs
			#for each iteration
			for(.iter in 1:dim(out_trim$Zs)[2]){
				
				#Store current states
				Znow<-factor(out_trim$Zs[,.iter])    
				
				#identify potential candidate switches:
				CandiCells<-table(Znow , Zref)/apply(table(Znow , Zref), 1, sum)>PropMin
				getCandi<-function(x) { as.numeric(as.character(row.names(as.matrix(x))[x])) }
				ListCandi<- apply(CandiCells, 1, getCandi)
				
				#IF its a matrix, turn into list. If its a single choice, go fast.
				# R stuff to make sure it deals with inputs correctly
				if(class(ListCandi)=='matrix'){
				ListCandi<-split(ListCandi, rep(1:ncol(ListCandi), each = nrow(ListCandi)))
				Candies<-expand.grid(ListCandi)  # each row is a labelling
				names(Candies)<-row.names(CandiCells)   # RAAAAH
				} else if (class(ListCandi)=='numeric'){
				Candies<-ListCandi
				} else {
				Candies<-expand.grid(ListCandi)  # each row is a labelling
				names(Candies)<-row.names(CandiCells)   # RAAAAH
				}

				namesCandies<-names(Candies)
				# Catch if no appropriate seperation available at all, check all permutations. Only occurs in really bad models with bad convergence.				
				if(class(Candies)=='data.frame'){
				if  ( max(sapply(apply(Candies, 1, unique), length))<length(row.names(CandiCells))){
					Candies<- permutations(K)}
					}else{
						if  (length(unique(Candies))<length(row.names(CandiCells))){
					Candies<- permutations(K)}
					} 

				MinusRefPars<-function(x) 	{ flp<- as.numeric(  row.names(CandiCells)[unlist(Candies[x,])])
						if(length(unique(flp))<length(flp)) { Inf
						} else {sum(abs( (refComp	-  c(out_trim$P[.iter,flp], out_trim$Mu[.iter,flp],out_trim$Sig[.iter,flp]))/refComp))	}}
											
				if( sum( apply(CandiCells, 1, sum)) >  dim(CandiCells)[1] ){
					BestOne<-which.min( sapply(1:dim(Candies)[1] , MinusRefPars))  # find the best perm out of options
					BestOne<-Candies[BestOne,]
					} else {BestOne<- Candies }   # chose this one if no comparing needed
				if(is.null(names(BestOne))) {names(BestOne)<-namesCandies}
	
				# Allocations
				Znew<-Znow; levels(Znew)<-as.numeric(BestOne)
				Zfixed[,.iter]<-as.numeric(as.character(Znew))
				# Parameters
				combinePars<-cbind(.iter,as.numeric(BestOne),  out_trim$Ps[.iter,as.numeric(names(BestOne))],out_trim$Mu[.iter,as.numeric(names(BestOne))], out_trim$Sig[.iter,as.numeric(names(BestOne))] )[order(as.numeric(BestOne), decreasing=FALSE),]
				AllPars<-rbind(AllPars, combinePars)
				}
			names(AllPars)<-c('Iteration', 'k', 'P', 'Mu','Sig')
			
			# sumarise Zs (find max)
			maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
			Zhat<- factor( apply(t(Zfixed), 2,maxZ))
			levels(Zhat)<- levels(Zhat)<-as.character(BestOne)
			return(list('Pars'=AllPars, 'Zs'=Zfixed))
			}
		