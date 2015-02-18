#' This functionZmix_multi_tempered
#' @param stuff and more
#' @keywords multi
#' @export
#' @examples
#' #not run
	Zmix_multi_tempered<-function(YZ, iter, k, alphas, sim=TRUE, EndSize=500){
					
					dMvn <- function(X,mu,Sigma) {
						    k <- ncol(X)
						    rooti <- backsolve(chol(Sigma),diag(k))
						    quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
						    return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))}
					trimit<-function(Out=Out, nEnd=EndSize){
							yo<-length(Out$Mu)		
							ps<-Out$P[[yo]][c(iter-nEnd+1):iter,]			
							mu<-subset(Out$Mu[[yo]], Iteration>c(iter-nEnd))
							covs<-subset(Out$Cov[[yo]], Iteration>c(iter-nEnd))
							zs<-Out$Zs[[yo]][c(iter-nEnd+1):iter,]
							Loglike<-Out$Loglike[c(iter-nEnd+1):iter]
			     			SteadyScore<-Out$SteadyScore$K[c(iter-nEnd+1):iter]
							list(Mu = mu,Cov=covs, P= ps,  Zs=zs, Y=Out$Y, Loglike=Loglike, SteadyScore=SteadyScore)	}
					dDirichlet<-function (x, alpha, log = FALSE) {
						    dlog = lgamma(sum(alpha)) + sum((alpha - 1) * log(x)) - sum(lgamma(alpha))
						    result = ifelse(!log, exp(dlog), dlog)
						    return(result)
							}

					parallelAccept<-function(w1, w2, a1, a2){
							w1[w1< 1e-200]<-1e-200   # truncate so super small values dont crash everyting
							w2[w2< 1e-200]<-1e-200
								T1<-dDirichlet(w2, a1, log=TRUE)
								T2<-dDirichlet(w1, a2, log=TRUE)
								B1<-dDirichlet(w1, a1, log=TRUE)
								B2<-dDirichlet(w2, a2, log=TRUE)
								MH<-min(1,	exp(T1+T2-B1-B2)) 
							Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
							return(Ax)}
					CovSample<-function(nsk, WkZk, ybark){
			 			ck<-c0+nsk/2
							if (nsk==0) {
								solve(rwish(c0, solve(C0)))
							} else {
							Ck<-	C0	+0.5*((nsk*n0)/(nsk+n0)*crossprod(ybark-b0)) +0.5*WkZk
							solve(rwish(ck, solve(Ck)))}}
					MuSample<-function(newCovLISTk, nsk,WkZk, ybark){
							newCovLISTk<-matrix(newCovLISTk, nrow=r, byrow=T)

							  if (nsk==0) {
							 rmvnorm(1, b0, newCovLISTk/n0)
							  } else {
								bk<-(n0/(nsk+n0))*b0+(nsk/(n0+nsk))*ybark
								Bk<-(1/(nsk+n0))*newCovLISTk    
							rmvnorm(1, t(bk), Bk)

								}}		
					minions<-function(ZZ){  # intake classifications
						# ns	
						IndiZ <- (ZZ == matrix((1:k), nrow = n, ncol = k, byrow = T))			
						ns <- apply(IndiZ,2,sum)  #  size of each group
							
							.Ysplit<-replicate(k, list())	#storage to group Y's
							WkZ<-	 replicate(k, list(0))	#storage for within group variability
							ybar<-	 replicate(k, list(0))
							for (.i in 1:k){					# grouping y's by Zs
								.Ysplit[[.i]]<-Y[ZZ==.i,]
								if (ns[.i]>1){					# for groups with >1 obsevations
							ybar[[.i]]<- as.matrix(t(apply(.Ysplit[[.i]], 2, mean) ))
								} else if (ns[.i]==1){			# if n=1 mean =y 
							ybar[[.i]]<-	t( as.matrix(.Ysplit[[.i]]))
								} else {  ybar[[.i]]<-	NA}
							 
							 #Within group unexplained variability
								if (ns[.i]==0) { WkZ[[.i]]<-	NA	} else if (ns[.i]==1){
								WkZ[[.i]]<-crossprod(as.matrix(.Ysplit[[.i]]-ybar[[.i]]))
								} else {	for ( .n in 1:ns[.i]){ 
								WkZ[[.i]]<-WkZ[[.i]]+ crossprod( .Ysplit[[.i]][.n,]-ybar[[.i]])
							}}}

							list('ns'=ns,'ybar'=ybar, 'WkZ'=WkZ)}

					if(sim==TRUE){ Y<- as.matrix(YZ$Y) }else{Y<- as.matrix(YZ)}  	 # change as needed for sim VS case studies, could add an option in function
				
					nCh<-length(alphas)	
					r<-dim(Y)[2] ;	n <-dim(Y)[1] 			
					n0=1   # this is tau, right?
			      	Mus<- 	 replicate(k, list())  ##
			      	Covs<- 	replicate(k, list()) ##   # k vector per iteration  (fix so as to save iterations)		
			  		 
					# storing final:
							Ps <-   replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
							Zs<- 	replicate(nCh,  matrix(0,nrow = iter, ncol = n)	, simplify=F)
							v<-data.frame(matrix(NA, nrow=0, ncol=r*r+2))
							FINcov<-replicate(nCh, v	, simplify=F) 
							v2<-data.frame(matrix(NA, nrow=0, ncol=r+2))
							FINmu<-replicate(nCh, v2	, simplify=F) 
							Loglike <-   matrix(0,nrow = iter, ncol = 1)
							SteadyScore<-data.frame("Iteration"=c(1:iter), "K0"=0)
					#hyperpars 	
					Ck<-	replicate(k, list())
					b0<-apply(Y,2,mean)
					c0<-2.5+(r-1)/2  ; if (c0<r) c0<-r+1		
					C0<-0.75*cov(Y)	
					d<-sum(c(1:r))+r
				
					# STEP 1: initiate groups   (iteration 1 only)

					pb <- txtProgressBar(min = 0, max = iter, style = 3)

					for (.it in 1:iter){  #for each iteration	
					  if(.it %% 10 == 0) { Sys.sleep(0.01)
					  	par(mfrow=c(2,1))
	  					plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
	  					ts.plot(Ps[[nCh]], main='emptying', col=rainbow(k))
					  Sys.sleep(0)
						setTxtProgressBar(pb, .it) }

					for (.ch in 1:nCh){   #for each chain

						if (.it==1){ 
							.Zs <-as.vector(kmeans(Y, centers=k)$cluster)
							bits<-minions(.Zs)
							ns<-bits$ns; ybar<-bits$ybar 
							WkZ<-bits$WkZ

							ns <-apply((.Zs  ==  matrix((1:k), nrow = n, ncol = k, byrow = T)),2,sum)
							} else {  # uses the allocations from end of last iteration
							bits<-minions(Zs[[.ch]][.it-1,])
								ns<-bits$ns; ybar<-bits$ybar ;	WkZ<-bits$WkZ
								.Zs<-Zs[[.ch]][.it-1,]
							}
				
				
						# STEP 2.1 : GENERATE Samples for WEIGHTS from DIRICHLET dist
								#  COUNTER FOR ALPHA CHANGE   if (ShrinkAlpha==TRUE){	
						
								Ps[[.ch]][.it,] = rdirichlet(m=1,par= ns+alphas[.ch])
							
						#STEP 2.2 GENERATE Samples from Covariance Matrix for each component
							newCov<-mapply( CovSample, ns, WkZ, ybar)
							newCovLIST<-as.list(as.data.frame(newCov))

							 FINcov[[.ch]]<-rbind(FINcov[[.ch]],cbind(t(newCov), 'K'=1:k, 'Iteration'=.it))


						#STEP 2.3 GENERATE SAMPLEs of the component specific Mu's from multivariate normal (bk,Bk)
								newMu<-mapply(MuSample, newCovLIST, ns, WkZ , ybar)
								newMuLIST<-as.list(as.data.frame(newMu))
								FINmu[[.ch]]<-rbind(FINmu[[.ch]],cbind(t(newMu), 'K'=1:k, 'Iteration'=.it))


						# STEP 3: Draw new classification probabilities:
									
									PZs<-matrix(0,ncol=k, nrow=n)	
								for (i in 1:n) { for (.k in 1:k){
								PZs[,.k]<-
								#dmvnorm(Y, newMuLIST[[.k]], matrix(newCovLIST[[.k]], nrow=r, byrow=T))*Ps[[.ch]][.it,.k]
								dMvn(Y, newMuLIST[[.k]], matrix(newCovLIST[[.k]], nrow=r, byrow=T))*Ps[[.ch]][.it,.k]
								}}

									#scale each row to one
									for (i in 1:n)		{
										if (sum(PZs[i,])==0) {
											PZs[i,]<-rep(1/k,k)    # if all probs are zero, randomly allocate obs. very rare, might lead to crap results
											}else {
									PZs[i,]<-PZs[i,]/sum(PZs[i,])}}
									## NEW allocations based on probabilities 
								for (i in 1:n){	Zs[[.ch]][.it,i]=sample((1:k),1, prob=PZs[i,])}	
							
							} # end of chain loop
					SteadyScore$K0[.it]<-sum(table(Zs[[nCh]][.it,])>0)

					## PARALLEL TEMPERING MOVES 
					 if(.it>20){
					 if( sample(c(1,0),1, 0.5)==1){	 
					 	Chain1<-sample( 1:(nCh-1), 1) ; 	Chain2<-Chain1+1
						MHratio<- parallelAccept(Ps[[Chain1]][.it,], Ps[[Chain2]][.it,], rep(alphas[Chain1],k), rep(alphas[Chain2],k))
						if (MHratio==1){  
						# Just flip the allocations since all pars drawn from this 
							.z1<-	Zs[[Chain1]][.it,] 	;		.z2<-	Zs[[Chain2]][.it,]
							Zs[[Chain1]][.it,]<-.z2		;		Zs[[Chain2]][.it,]<-.z1
							}}}
					
					                                                                                                             #logLikelihood
						for (i in 1:n){
						non0id<-c(1:k)[ns > 0]
							.ll<-0
							for (numK in 1:length(non0id)){
								 .ll<-.ll+ Ps[[nCh]][.it,non0id[numK]]* dmvnorm(Y[i,], newMu[,non0id[numK]], matrix(newCov[,non0id[numK]], ncol=r,nrow=r, byrow=T))
							}
						Loglike[.it]<-Loglike[.it]+ log(.ll)}	
						
					}
					close(pb)

					 # end of iterations loop

					bigres<-list( P=Ps, Cov=FINcov, Mu=FINmu, Zs=Zs, YZ=YZ, SteadyScore=SteadyScore, Loglike=Loglike)
					return(trimit(bigres, nEnd=EndSize))
					 }