#' Label switch fix for Multivariate Normals 
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #not run

	QuickSwitchMVN<-function(MVNResult,PropMin=0.3){
			K<-dim(MVNResult$P)[2];r<-dim(MVNResult$Y$Y)[2]
			wml<-which.max(MVNResult$Loglike)
			Zref<-factor(MVNResult$Zs[wml,], levels=1:K,ordered=FALSE) 
			# if (LineUpBy==3){
			# 			FinalOrderChoice<-order(MVNResult$P[wml,], decreasing=TRUE)
			# } else if(LineUpBy==4){
			# 			FinalOrderChoice<-order(MVNResult$P[wml,], decreasing=TRUE)}
			FinalOrderChoice<-order(MVNResult$P[wml+1,], decreasing=TRUE)  # gets too complex to use other pars so limiting to weights for now		
			levels(Zref)<- c(1:K)[order(FinalOrderChoice)]
			Zref<- factor(Zref,levels(Zref)[order(levels(Zref))])
			
			## UPDATE HERE 
			non0ref<-sum(table(Zref)>0)			##
			iterationID<- max(MVNResult$Mu$Iteration)-(dim(MVNResult$P)[1] - wml)
			##REF PARAMETERS FOR SECOND PHASE
			P_ref<- MVNResult$P[wml+1,FinalOrderChoice[1:non0ref]]		
			Mu_ref<-MVNResult$Mu[MVNResult$Mu$Iteration==iterationID+1,]
			Mu_ref<-Mu_ref[FinalOrderChoice[1:non0ref],]
			Mu_ref<-  stack( Mu_ref[,-c(dim(Mu_ref)[2]-1, dim(Mu_ref)[2])])[,1]

			refPar<-c(P_ref, Mu_ref)	##

			Zfixed<-MVNResult$Zs[c(1:dim(MVNResult$Zs)[1]-1),]
			iterID1<-MVNResult$Mu$Iteration[1]+1
			AllPars<-data.frame(diag(3+r+(r*r)))[numeric(0), ]

			for(.iter in 1:(dim(MVNResult$Zs)[1]-1)){   ###
							Znow<-factor(MVNResult$Zs[.iter,])  
							CandiCells<-table(Znow , Zref)/apply(table(Znow , Zref), 1, sum)>PropMin
							getCandi<-function(x) { as.numeric(as.character(row.names(as.matrix(x))[x])) }
							ListCandi<- apply(CandiCells, 1, getCandi)


			#IF its a matrix, turn into list. If its a single choice, go fast
				if(class(ListCandi)=='numeric'){			Candies<-t(as.matrix(ListCandi, byrow=T))
				} else if (class(ListCandi)=='matrix'){		ListCandi<-split(ListCandi, rep(1:ncol(ListCandi), each = nrow(ListCandi)))
															Candies<-expand.grid(ListCandi)  
															names(Candies)<-row.names(CandiCells)   
				} else {									Candies<-expand.grid(ListCandi)  # each row is a labelling
															names(Candies)<-row.names(CandiCells)   # RAAAAH
				}
				
				
				if(dim(Candies)[1] >  1){Candies<-Candies[apply(Candies, 1, function(x) sum(table(x)>0))==dim(Candies)[2],]}

				#if( sum( apply(CandiCells, 1, sum)) >  dim(CandiCells)[1]){ 
				
					#refPar<-refPar[1:dim(Candies)[2]]}  # store reference pars only for not emoty groups
				MinusRefPars<-function(x) { 
						flp<- as.numeric(  names(Candies)[unlist(Candies[x,])]); 
					if(length(unique(flp))<length(flp)){ "NA"
					} else{
						comparePar_p<-MVNResult$P[.iter+1,flp] 
						CurrentiterationID<- max(MVNResult$Mu$Iteration)-(dim(MVNResult$P)[1]) +(.iter+1)				
						comparePar_Mu<-MVNResult$Mu[MVNResult$Mu$Iteration==CurrentiterationID,]
						comparePar_Mu<-comparePar_Mu[flp,]
						comparePar_Mu<-  stack( comparePar_Mu[,-c(dim(comparePar_Mu)[2]-1, dim(comparePar_Mu)[2])])[,1]
						currentPars<-c(comparePar_p, comparePar_Mu)
						sum(abs((refPar-currentPars)/refPar ))}}

					
					if(dim(Candies)[1] >  1){
							BestOne<-which.min( sapply(1:dim(Candies)[1] , MinusRefPars))  # find the best perm out of options
							BestOne<-Candies[BestOne,]
						} else {BestOne<- Candies }   # chose this one if no comparing needed

					#now move everything according to choice
					Znew<-Znow; levels(Znew)<-as.numeric(BestOne)
					Zfixed[.iter,]<-as.numeric(as.character(Znew))
					# Parameters

					swP<-cbind(.iter,1:K,  MVNResult$P[.iter+1,])
					swP<-swP[as.numeric(colnames(BestOne)),]
					colnames(swP)<-c('Iteration', 'k', 'P')
					#MU
					swM<-MVNResult$Mu[MVNResult$Mu$Iteration==iterID1,]
					swM<-swM[as.numeric(colnames(BestOne)),1:r]
			        colnames(swM)<-paste("Mu", 1:r, sep='_')
			        swCV<-MVNResult$Cov[MVNResult$Cov$Iteration==iterID1,]
					swCV<-swCV[as.numeric(colnames(BestOne)),1:(r*r)]
					colnames(swCV)<-paste("Cov", 1:(r*r), sep='_')


					combinePars<- cbind(swP, swM, swCV)
					iterID1<-iterID1+1
					#combinePars<-combinePars[as.numeric(names(BestOne)),]
					combinePars[,2]<-as.numeric(BestOne)

					AllPars<-rbind(AllPars, combinePars)
				}
			 		
						maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
						Zhat<- factor( apply(t(Zfixed), 1,maxZ), levels=1:K)
								
						varSum<- sum(apply(AllPars[,-c(1,2)], 2, var)) 

				return(list(Pars=AllPars, Zs=Zfixed, YZ=MVNResult$Y, varSum=varSum, SteadyScore=MVNResult$SteadyScore))
			}
