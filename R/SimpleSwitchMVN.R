#' This function fixes the label switching in the output of a Multivariate Gibbs sampler. Requires a constant number of non-empty groups in samples, but can deal with empty groups.
#' @param pickRef;      'random'= reference picked as a random iteration
#' 			'last'= Reference is final iteration
#' 			"MaxLike"= Reference is chosen as iteration maximising Log LIkelihood,  (MVNResult$Loglike)
#' @param 	Numdim: NMumber of dimensions
#' @param 	MVNResult: LIST of output of Gibbs sampler with the following components:
#' 
#' 		[[Zs]] estimated component (ALLOCATIONS) for each Y at each iteration (by rows=Iterations, Columns=1:n)  
#'
#' 		[P]] matrix of WEIGHTS, columns - K, rows = Iterations
#'
#' 		[[Mu]] Data frame of MEANS  with columns: Mu_1. Mu_2, .. ...,Mu_dim ,  k (this is the group), Iteration.  
#'
#' 		[[Loglike]] OPTIONAL: Loglikelihood (vector the length of iterations). Needed only if PickRef=="MaxLike"
#'

#' @keywords Wishart
#' @export
#' @examples
#' RESULT$
#' 	Unswitched<-SimpleSwitchMVN( RESULT, pickRef="last")

	SimpleSwitchMVN<-function(MVNResult, pickRef=c("random", "last", "MaxLike"), NumDim=2, PropMin=0.01){
		permutations <- function(n){
			if(n==1){
			return(matrix(1))
			} else {
			sp <- permutations(n-1)
			p <- nrow(sp)
			A <- matrix(nrow=n*p,ncol=n)
			for(i in 1:n){
			A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
			}
			return(A)
			}}
		getCandi<-function(x) { as.numeric(as.character(row.names(as.matrix(x))[x])) }
		maxZ<-function (x)  as.numeric(names(which.max(table( x ))))

		K<-dim(MVNResult$P)[2]
		rdim<-NumDim		# number of dims
 
		# pick 
		if(pickRef=="random"){wml<- sample(1:dim(MVNResult$Zs)[1]) 
		} else if (pickRef=="last") {  wml<-dim(MVNResult$Zs)[1]
		} else if (pickRef--"MaxLike"){wml<-which.max(MVNResult$Loglike)}

		Zref<-factor(MVNResult$Zs[wml,], levels=1:K,ordered=FALSE) 
		non0ref<-sum(table(Zref)>0)
		
		FinalOrderChoice<-order(MVNResult$P[wml,], decreasing=TRUE)[1:non0ref]  #USE weights

		Zref<-factor(Zref) 
		levels(Zref)<- c(1:K)[order(FinalOrderChoice)]
		Zref<- factor(Zref, levels=levels(Zref)[order(levels(Zref))])  # tidy up factor levels
# rename iterations
 MVNResult$Mu$Iteration<-as.vector(sapply( c(1: (length(MVNResult$Mu$Iteration)/K) ),function(x) rep(x, K))	)
# MVNResult$Cov$Iteration<-as.vector(sapply( c(1: (length(MVNResult$Cov$Iteration)/K) ),function(x) rep(x, K))	)
			# REF PARAMETERS FOR SECOND PHASE
			P_ref<- MVNResult$P[wml,FinalOrderChoice]		
			Mu_ref<-MVNResult$Mu[MVNResult$Mu$Iteration==wml,1]
			Mu_ref<-Mu_ref[FinalOrderChoice]
			refPar<-c(P_ref, Mu_ref)	##

			Zfixed<-matrix(data=NA, nrow=dim(MVNResult$Zs)[1], ncol=dim(MVNResult$Zs)[2])
			AllPars<-data.frame(diag(3+rdim+(rdim*rdim)))[numeric(0), ]

		for(.iter in 1:(dim(MVNResult$Zs)[1])){   
				
#iterIDnow<-iterID1-1+.iter
			Znow<-factor(MVNResult$Zs[.iter,])  
			CandiCells<-table(Znow , Zref)/apply(table(Znow , Zref), 1, sum)>PropMin
			ListCandi<- apply(CandiCells, 1, getCandi)

		#IF its a matrix, turn into list. If its a single choice, go fast
			if(class(ListCandi)=='numeric'){
				Candies<-t(as.matrix(ListCandi, byrow=T))
			} else if (class(ListCandi)=='matrix'){		
				ListCandi<-split(ListCandi, rep(1:ncol(ListCandi), each = nrow(ListCandi)))
				names(ListCandi)<-row.names(CandiCells)   
				Candies<-expand.grid(ListCandi)  
			} else {	
				Candies<-as.data.frame(expand.grid(ListCandi))  # each row is a labelling
				names(Candies)<-row.names(CandiCells)   
			}
			
			if(dim(Candies)[1] >  1){
				NumberValid<-sum(apply(Candies, 1, function(x) sum(table(x)>0))==dim(Candies)[2])
				if (NumberValid>0){
					Candies<-Candies[apply(Candies, 1, function(x) sum(table(x)>0))==dim(Candies)[2],]
				} else{
				Candies<-matrix(as.numeric(colnames(CandiCells)[t(permutations(non0ref))]), ncol=non0ref, byrow='TRUE')
				colnames(Candies)<- rownames(CandiCells)
					}
}



if(dim(Candies)[1] >  1){


countDiffs<-rep(Inf, dim(Candies)[1])

	now<-na.omit( as.numeric(  row.names(CandiCells)))
	cmu<-MVNResult$Mu[MVNResult$Mu$Iteration==.iter,1]
	
for (i in 1:dim(Candies)[1] ){
	preswap<-data.frame( "K"=as.numeric(  row.names(CandiCells)),  MVNResult$P[.iter,now], cmu[now])
# replace K column with proposed values, reorder and turn into vector.
preswap$K<-as.numeric(Candies[i,])
preswap<-preswap[ order(preswap$K, decreasing=FALSE),]
currentPars<-stack(preswap[,-1])[,1]
countDiffs[i]<-sum(abs((refPar-currentPars)/refPar ))
}					
BestOne<-Candies[which.min(countDiffs),]

} else {	BestOne<- as.data.frame( Candies) }   # chose this one if no comparing needed

	#now move everything according to choice
					Znew<-Znow
					levels(Znew)<-as.numeric(BestOne)
					Zfixed[.iter,]<-as.numeric(as.character(Znew))
					# Parameters
				swP<-cbind(.iter,1:K,  MVNResult$P[.iter,])
				colnames(swP)<-c('Iteration', 'K', 'P')

				swM<-MVNResult$Mu[MVNResult$Mu$Iteration==.iter, -(rdim+2) ] 
				colnames(swM)<-c(paste("Mu", 1:rdim, sep='_'), "K")
				#swCV<-MVNResult$Cov[MVNResult$Cov$Iteration==.iter, -(rdim*rdim+2) ]   
				#colnames(swCV)<-c(paste("Cov", 1:(rdim*rdim), sep='_'), "K")
				#combinePars<-merge(merge(swP, swM, by="K"), swCV, by="K")	#extract non-empty
				combinePars<-merge(swP, swM, by="K")	#extract non-empty
				combinePars<-combinePars[combinePars$K %in% as.numeric(names(BestOne)) ,]
				#rename to match labels
				combinePars$K<-as.numeric(BestOne)
				
				AllPars<-rbind(AllPars, combinePars[order(combinePars$K),])
				}
	
					Zhat<- factor( apply(t(Zfixed), 1,maxZ), levels=1:K)
								
					varSum<- sum(apply(AllPars[,-c(1,2)], 2, var)) 

				return(list(Pars=AllPars, Zs=Zfixed))
			}
