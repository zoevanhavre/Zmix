
#' Run Gibbs sampler with prior tempering for univariate gaussian mixture models
#'
#' This function ...
#' @param y, k,iter, isSim=FALSE, alphas
#' @keywords Gibbs sampler, univariate, normal, gaussian, mixture, prior parallel tempering.
#' @export
#' @examples
#' # Coming zoon

Zmix_univ_tempered<-function(y, k=10,iter=5000,  tau=0.01, isSim=FALSE, alphas= c(30, 20, 10, 5, 3, 1, 0.5, 1/2^(c(2,3,4,5,6, 8, 10, 15, 20, 30)))){
				if(isSim==TRUE) {Y<-y$Y
					}else{ Y<-y}
				parallelAccept<-function(w1, w2, a1, a2){
						w1[w1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
						w2[w2< 1e-200]<-1e-200
						T1<-dDirichlet(w2, a1, log=TRUE)
						T2<-dDirichlet(w1, a2, log=TRUE)
						B1<-dDirichlet(w1, a1, log=TRUE)
						B2<-dDirichlet(w2, a2, log=TRUE)
						MH<-min(1,	exp(T1+T2-B1-B2))
						Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
						return(Ax)}


					nCh<-length(alphas)
						TrackParallelTemp<-matrix(nrow=iter, ncol=nCh)
						TrackParallelTemp[1,]<-c(1:nCh)


					#tau=1
					n <-length(Y)
					a=2.5; b=2/var(Y)
					d<-2
					lambda=sum(Y)/n
			                                                                                                          # 1. set up priors
					mux<-list(mu=seq(from=min(Y), to=max(Y),length.out=k),sigma=rep(1, k),p=rep(1/k,k), k=k)
					n <-length(Y)
					a=2.5; b<-0.5*var(Y);d<-2
					lambda<-sum(Y)/n  ;
					pb <- txtProgressBar(min = 0, max = iter, style = 3)

			                                                                                                          # 2. set up matrices for parameters which will be saved
					map =    matrix(0,nrow = iter, ncol = 1)
					Loglike =   matrix(0,nrow = iter, ncol = 1)
					Bigmu = replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
					Bigsigma=replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
					Bigp =  replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
					Pzs =   replicate(nCh,  matrix(0,nrow = n, ncol = k)	, simplify=F)
					ZSaved=	replicate(nCh,  matrix(0,nrow = n, ncol = iter)	, simplify=F)
					SteadyScore<-data.frame("Iteration"=c(1:iter), "K0"=k)


			                                                                                                          # start chains and  create inits needed for i=1
					for (.ch in 1:nCh){
					Bigmu[[.ch]][1,] <- mux$mu                                                                        # initial value of mu's
					mu0=mux$mu
					Bigp[[.ch]][1,] = mux$p                                                                           # inital value of p's
					p0=mux$p
					Bigsigma[[.ch]][1,] = mux$sigma                                                                   # inits for sigma
					sig0=mux$sigma	}

			                                                                                                          #  		initialize chains:  ie. iteration 1:
					j<-1

			                                                                                                          # 1. Generate P(Z[i](t)=j)  for j=1,...,known
			                                                                                                          # from paper, computing pzs
					for (.ch in 1:nCh){
					for (i in 1:n) {
					Pzs[[.ch]][i,]<-(p0/sqrt(sig0))*exp(-((Y[i]-mu0)^2)/(2*sig0))

					Pzs[[.ch]][i,]<-Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,]) }}											  # Scale to equal 1?

			                                                                                                          # 2 Make indicator matrix of assignments based on Pzs
			                                                                                                          #	sample 1 of the k classes for each row by Pzs (prob)
					for (.ch in 1:nCh){
					Z<-matrix()
					for (i in 1:n){Z[i]=sample((1:k),1, prob=Pzs[[.ch]][i,])}
					matk = matrix((1:k), nrow = n, ncol = k, byrow = T)
					IndiZ = (Z == matk)
					ZSaved[[.ch]][,1]<-Z
			                                                                                                          # 3 compute ns and sx
					ns = apply(IndiZ,2,sum)
					for (i in 1:length(ns)){ if ( is.na(ns[i])) ns[i]<-0 }
					sx = apply(IndiZ*Y, 2, sum)

			                                                                                                          # 4 Generate P[j](t) from dirichlet  (and save)
					Bigp[[.ch]][j,] = rdirichlet(m=1,par= ns+alphas[.ch])

			                                                                                                          # 5	Generate Mu's   (and save)
					Bigmu[[.ch]][j,]<-rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][1,]/(tau+ns))) # must be sqrt as r takes in sd not var

			                                                                                                          #	ybar<-sx/ns
					for (i in 1:length(Bigmu[[.ch]][j,])){ if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }
			                                                                                                          # 6  Compute sv[j](t)
					.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
					for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
					sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)                                                         # changes, added /ns

			                                                                                                          # 7 Generate Sigma's (and save)

					Bigsigma[[.ch]][j,]<- rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)

					}

			                                                                                                          # Log Likelihood: # Sum-n (log Sum-K ( weights x dnorm (y,thetas)))
					for (i in 1:n){
					non0id<-c(1:k)[ns > 0]
					Loglike[j]<-Loglike[j]+ log(
						 sum( Bigp[[nCh]][j,non0id]*dnorm(Y[i], mean=Bigmu[[nCh]][j,non0id], sd=sqrt(Bigsigma[[nCh]][j,non0id]))))}

			                                                                                                          ### now finish loop for j>2


					for (j in 2:iter){


					 if(j %% 100==0){
					Sys.sleep(0.01)
					setTxtProgressBar(pb, j)
					# par(mfrow=c(1,3))
					# plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
					# ts.plot(Bigp[[nCh]], main='Weights from target posterior', col=rainbow(k))
					# ts.plot(TrackParallelTemp[,c(nCh:1)], main='Track Parallel Tempering', col=rainbow(nCh))
					 	Sys.sleep(0)}

					for (.ch in 1:nCh){
			                                                                                                          #################

			                                                                                                          # 1. Generate P(Z[i](t)=j)
					for (i in 1:n) {
					Pzs[[.ch]][i,]<-(Bigp[[.ch]][j-1,]/sqrt(Bigsigma[[.ch]][j-1,]))*exp(-((Y[i]-Bigmu[[.ch]][j-1,])^2)/(2*Bigsigma[[.ch]][j-1,]))
			                                                                                                          # Scale to equal 1
					Pzs[[.ch]][i,]<-Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,])
					}
			                                                                                                          # 2 Make indicator matrix of assignments based on Pzs
			                                                                                                          #	sample 1 of the k classes for each row by Pzs (prob)
					for (i in 1:n){Z[i]=sample((1:k),1,replace=T, prob=Pzs[[.ch]][i,])}
					matk = matrix((1:k), nrow = n, ncol = k, byrow = T)
			                                                                                                          #indicator
					IndiZ = (Z == matk)
					ZSaved[[.ch]][,j]<-Z
			                                                                                                          # 3 # compute ns and sx
					ns = apply(IndiZ,2,sum)
			                                                                                                          # fix if ns=NA to =0
					for (i in 1:length(ns)){
					if ( is.na(ns[i])) ns[i]<-0 }

					sx = apply(IndiZ*Y, 2, sum)

			                                                                                                          # 4 Generate P[j](t) from dirichlet  (and save)
					Bigp[[.ch]][j,] = rdirichlet(m=1,par= ns+alphas[.ch])

			                                                                                                          # 5	Generate Mu's   (and save)
					Bigmu[[.ch]][j,]<-rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][j-1,]/(tau+ns)))
					for (i in 1:length(Bigmu[[.ch]][j,])){ if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }

			                                                                                                          # 6 Compute sv[j](t)
					.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
					for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
					sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)                                                         # changes, added /ns

			                                                                                                          # 7 Generate Sigma's (and save)
					Bigsigma[[.ch]][j,]<-  rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)

			                                                                                                          #calculate MAP estimates at this iteration:
			                                                                                                          #if (.ch==nCh){
			                                                                                                          #mapy[j,]<- prod(apply(Pzs[[.ch]], 2,sum))
			                                                                                                          #mapmu<-prod(dnorm(Bigmu[[.ch]][j,], mean=lambda, sd=(Bigsigma[[.ch]][j,]/tau) ))
			                                                                                                          #mapsig<- prod(dinvgamma(Bigsigma[[.ch]][j,], a,b))
			                                                                                                          #map[j,]<-prod(mapy[j,], mapmu, mapsig)}
					}
			           #new                                                                                               #### PARALLEL TEMPERING MOVES ###

					 if(j>1 ) {TrackParallelTemp[j,]<-TrackParallelTemp[j-1,]}      # SET para chains to previous values                                                                                           # how often??? lets go with probability 10% of switching at any time

					if(j>20 && nCh>1){
					if( sample(c(1,0),1, 0.9)==1){		# FREQ OF TEMPERING!
			                                                                                                          #Pick chains, just chose one and the one next to it
      if( j%%2==0){chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2==0]   #evens
      } else {chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2!=0]  }   #odds
	if(nCh==2){ chainset<-c(1) }
 for( eachChain in 1:length(chainset)){
                Chain1<-chainset[eachChain]
                Chain2<-Chain1+1

					#Chain1<-sample( 1:(nCh-1), 1)
					#Chain2<-Chain1+1
			       ## allow non-adjacent chains
					#Chain1<-sample( c(1:nCh), 1)
					#Chain2<-sample(c(1:nCh)[c(1:nCh)!=Chain1],1)
			                                                                                                         # check ratio
					MHratio<- parallelAccept(Bigp[[Chain1]][j,], Bigp[[Chain2]][j,], rep(alphas[Chain1],k), rep(alphas[Chain2],k))
					if (MHratio==1){                                                                                  # switch the chains (weights, mean, sigma, Zs)
			             #new
			             .tpt1<-  TrackParallelTemp[j,Chain1 ]
			             .tpt2<-  TrackParallelTemp[j,Chain2 ]
						TrackParallelTemp[j,Chain1 ]<-.tpt2
			            TrackParallelTemp[j,Chain2 ]<-.tpt1
			                                                                                                   # Weights
					.p1<-	Bigp[[Chain1]][j,]
					.p2<-	Bigp[[Chain2]][j,]
					Bigp[[Chain1]][j,]<-.p2
					Bigp[[Chain2]][j,]<-.p1
			                                                                                                          # Means
					.m1<-	Bigmu[[Chain1]][j,]
					.m2<-	Bigmu[[Chain2]][j,]
					Bigmu[[Chain1]][j,]<-.m2
					Bigmu[[Chain2]][j,]<-.m1
			                                                                                                          # SD
					.s1<-	Bigsigma[[Chain1]][j,]
					.s2<-	Bigsigma[[Chain2]][j,]
					Bigsigma[[Chain1]][j,]<-.s2
					Bigsigma[[Chain2]][j,]<-.s1

			                                                                                                          # Zs
					.z1<-	ZSaved[[Chain1]][,j]
					.z2<-	ZSaved[[Chain2]][,j]
					ZSaved[[Chain1]][,j]<-.z2
					ZSaved[[Chain2]][,j]<-.z1
					}		}		}
					}
			                                                                                                          #logLikelihood
					for (i in 1:n){
					non0id<-c(1:k)[ns > 0]
					Loglike[j]<-Loglike[j]+ log( sum( Bigp[[nCh]][j,non0id]*dnorm(Y[i], mean=Bigmu[[nCh]][j,non0id], sd=sqrt(Bigsigma[[nCh]][j,non0id]))))}

					SteadyScore$K0[j]<-sum(table(ZSaved[[nCh]][,j])>0)

					}

					close(pb)

					# make at least one plot showing target
					Sys.sleep(1)
					par(mfrow=c(1,3))
					plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
					ts.plot(Bigp[[nCh]], main='Weights from target posterior', col=rainbow(k))
					ts.plot(TrackParallelTemp[,c(nCh:1)], main='Track Parallel Tempering', col=rainbow(nCh))

					#SteadyScore<-unlist(lapply(apply(ZSaved[[nCh]], 2, table), length))
					BigRes<-list(Bigmu = Bigmu, Bigsigma=Bigsigma, Bigp = Bigp, Loglike=Loglike, Zs=ZSaved, YZ=y, SteadyScore=SteadyScore,TrackParallelTemp=TrackParallelTemp)
					#SmallRes<-trimit(BigRes,  nEnd=EndSize)
					#return(SmallRes)
					return(BigRes)
					}
