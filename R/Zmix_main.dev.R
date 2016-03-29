library("mvtnorm")
library("dplyr")
library("compiler")
library("microbenchmark")
library("ggplot2")
library("mvnfast")

Zmix_main<-function(
          y,
          iterations=2000,
          k=10,
          alphas= c(1/2^(c(6, 10, 30))),
          burn=500,
          init.method="Kmeans",
          verbose=TRUE
          ){

      ## INNER FUNCTIONS
      # Initial clustering
      initiate_Z<-function(x, method=c("Kmeans", "random", "single") ){
        if(method=="Kmeans"){
        init.allocations<-  as.vector(kmeans(x, centers=k)$cluster)
        } else if (method=="random"){
          init.allocations<-  base::sample(c(1:k), n, replace=TRUE)
        } else if (method=="single"){
          init.allocations<-  rep(1, n)
        }
        return(init.allocations)
      }
      initiate_values<-function(ZZ){
        IndiZ <-  (ZZ == matrix((1:k), nrow = n, ncol = k, byrow = T))
        ns <-     apply(IndiZ,2,sum)  	#  size of each group
        .Ysplit<- replicate(k, list())	#storage to group Y's
        WkZ<-	    replicate(k, list(0))	#storage for within group variability
        ybar<-	  replicate(k, list(0))
        for (.i in 1:k){
        if(r>1){  .Ysplit[[.i]]<-y[ZZ==.i,]
        }else{ .Ysplit[[.i]]<-y[ZZ==.i]}

          if (ns[.i]>1){					# for groups with >1 obsevations
            if(r>1){
              ybar[[.i]]<- as.matrix(t(apply(.Ysplit[[.i]], 2, mean) ))  # UPDATE FOR UNIV
            } else {
            #  ybar[[.i]]<- as.matrix(t(apply(.Ysplit[[.i]], 2, mean) ))  # UPDATE FOR UNIV

            }
            } else if (ns[.i]==1){
              ybar[[.i]]<-	t( as.matrix(.Ysplit[[.i]]))
              } else {
                ybar[[.i]]<-	NA
              }
          #Within group unexplained variability
          if (ns[.i]==0) {
            WkZ[[.i]]<-	NA
            } else if (ns[.i]==1){
              WkZ[[.i]]<-crossprod(as.matrix(.Ysplit[[.i]]-ybar[[.i]]))
              } else {
                for ( .n in 1:ns[.i]){
                  WkZ[[.i]]<-WkZ[[.i]]+ crossprod( .Ysplit[[.i]][.n,]-ybar[[.i]])
                }
              }
            }
            list('ns'=ns,'ybar'=ybar, 'WkZ'=WkZ)
        }
      MuSample<-function(covs, ns_mu, WkZ_mu, ybar_mu){
          # covs<-matrix(covs, nrow=r, byrow=T)
          if (ns_mu==0) {
            rmvn(1, b0, covs/n0)
            } else {
              bk<-(n0/(ns_mu+n0))*b0+(ns_mu/(n0+ns_mu))*ybar_mu
              Bk<-(1/(ns_mu+n0))*covs
              rmvn(1, t(bk), Bk)
            }
        }
      Update_probs<-function(Mu, Cov, Pi ){
          update_group_Prob <-function(x){
              a1<-dmvn(y, Mu[[x]], Cov[[x]])
              a2<-a1*Pi[x]
              # for (i in 1:n)		{
              #   if (sum(PZs[i,])==0) {
              #     PZs[i,]<-rep(1/k,k)    # if all probs are zero, randomly allocate obs. very rare, might lead to crap results
              #   }}
              return(a2)
            }
          # apply to each component
          ugp<-mapply(update_group_Prob, c(1:k))
          # scale each row
          ugp/apply(ugp, 1, sum)
        }
      CovSample<-function(ns_cov, WkZ_cov, ybar_cov){
          c_cov<-c0+ns_cov
          if (ns_cov==0) {
            MCMCpack::riwish(c0, C0)
          } else {
            C_cov<- C0 +((ns_cov*n0)/(ns_cov+n0)*crossprod(ybar_cov-b0)) +WkZ_cov
            MCMCpack::riwish(c_cov,C_cov)
          }
        }
      Update_Zs<-function(PZs){mapply(function(x) sample( (1:k), 1, prob=PZs[x,]), c(1:n))}
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
      dDirichlet<-function (x, alpha, log = FALSE) {
          dlog = lgamma(sum(alpha)) + sum((alpha - 1) * log(x)) - sum(lgamma(alpha))
          result = ifelse(!log, exp(dlog), dlog)
          return(result)
      	}
      # function
 rdirichlet<-cmpfun(rdirichlet)

## compute basic values
  nCh<- length(alphas)
  if(is.vector(y)){
    # UNIVARIATE
    r<-1
    n<-length(y)
    ## univ hyper priors
    a<-       2.5
    b<-       2/var(y)
    lambda<-  sum(y)/n
    d<-   sum(c(1:r))+r
  } else {
    # MULTIVARIATE
    r<-   dim(y)[2]
    n<-   dim(y)[1]
    n0<-  1
    Ck<-	replicate(k, list())
    c0<-  r+1
    b0<-  apply(y,2,mean)
    C0<-  0.75*cov(y)
    d<-   sum(c(1:r))+r
  }

  Loglike <-   rep(0,iterations)
  SteadyScore<-data.frame("Iteration"=c(1:iterations), "K0"=0)
  map <-    matrix(0,nrow = iter, ncol = 1)




## parameters to estimate
## storage structure  par[[chain]][[iterations]][[k]]
  Ps<-     replicate(nCh, replicate(iterations, list()))
  Mus<-    replicate(nCh, replicate(iterations, list()))
  Covs<-   replicate(nCh, replicate(iterations, list()))
  Zs<-     replicate(nCh, replicate(iterations, list()))

  if(verbose==TRUE){pb <- txtProgressBar(min = 0, max = iterations, style = 3)}

# FOR EACH ITERATION
for (.it in 1:iterations){  #for each iteration
  # TRACKER
if(verbose==TRUE && .it %% 10 == 0) {
  Sys.sleep(0.01)
  if(r>1){
  par(mfrow=c(2,1))
  plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
  ts.plot( t(sapply(Ps[[nCh]], rbind)), main='Target Weights', col=rainbow(k))
  }
  Sys.sleep(0)
  setTxtProgressBar(pb, .it)
}

  for (.ch in 1:nCh){   #for each chain

# Input values
    if (.it==1){
      # iteration=1, initiallize input values
      init.Z<-  initiate_Z(y, init.method)
      input.values<-    initiate_values(init.Z) # compute values needed
        ns<-  input.values$ns
        ybar<-input.values$ybar
        WkZ<- input.values$WkZ
      } else {
      # Update given current pars
      input.values<-initiate_values(Zs[[.ch]][[.it-1]])  # check me
        ns<-    input.values$ns                         # check me
        ybar<-  input.values$ybar                       # check me
        WkZ<-   input.values$WkZ                        # check me
      }

      # STEP 2.1 : GENERATE Samples for WEIGHTS from DIRICHLET dist
      Ps[[.ch]][[.it]] <- rdirichlet(m=1, par= ns+alphas[.ch])

      # STEP 2.2 GENERATE Samples from Covariance
      Covs[[.ch]][[.it]]<-mapply(CovSample, ns, WkZ, ybar, SIMPLIFY=FALSE)

      # STEP 2.3 GENERATE SAMPLEs of Means
      Mus[[.ch]][[.it]]<- mapply(MuSample,Covs[[.ch]][[.it]], ns, WkZ, ybar, SIMPLIFY=FALSE)

      # STEP 3.1: Draw new classification probabilities:
      PZs<-Update_probs(Mus[[.ch]][[.it]], Covs[[.ch]][[.it]], Ps[[.ch]][[.it]])

      # STEP 3.2: Update allocations based on probabilitie
      Zs[[.ch]][[.it]]<-Update_Zs(PZs)

} # end chain loop


## PRIOR PARALLEL TEMPERING
if(.it>20 && runif(1)<.9){
  Chain1<-sample( 1:(nCh-1), 1) ; 	Chain2<-Chain1+1
  MHratio<- parallelAccept(Ps[[Chain1]][[.it]], Ps[[Chain2]][[.it]], rep(alphas[Chain1],k), rep(alphas[Chain2],k))
  if (MHratio==1){
  # Flip the allocations
  .z1<-	Zs[[Chain1]][[.it]] 	;.z2<-	Zs[[Chain2]][[.it]]
  Zs[[Chain1]][[.it]]<-.z2 	;Zs[[Chain2]][[.it]]<-.z1

  #Mu
  .mu1<-	Mus[[Chain1]][[.it]] 	;.mu2<-	Mus[[Chain2]][[.it]]
  Mus[[Chain1]][[.it]]<-.mu2 	;Mus[[Chain2]][[.it]]<-.mu1

  #Cov
  .cv1<-	Covs[[Chain1]][[.it]] 	;.cv2<-	Covs[[Chain2]][[.it]]
  Covs[[Chain1]][[.it]]<-.cv2 	;Covs[[Chain2]][[.it]]<-.cv1

  #Ps
  .p1<-	Ps[[Chain1]][[.it]] 	;.p2<-	Ps[[Chain2]][[.it]]
  Ps[[Chain1]][[.it]]<-.p2 	;Ps[[Chain2]][[.it]]<-.p1
}
}

# count number of occupied components
SteadyScore$K0[.it]<- factor(Zs[[.ch]][[.it]]) %>% levels(.) %>% length(.)

# compute Log Likelihood ( still need to optomize)
for (i in 1:n){
  non0id<-Zs[[nCh]][[.it]] %>% factor() %>% levels() %>% as.numeric()
    .ll<-0
    for (numK in 1:length(non0id)){
       .ll<-.ll+
       Ps[[nCh]][[.it]][non0id[numK]]*dmvn(
         y[i,], Mus[[nCh]][[.it]][[ non0id[numK] ]],
         Covs[[nCh]][[.it]][[ non0id[numK] ]])
    }
  Loglike[.it]<-Loglike[.it]+ log(.ll)
}


} # end iteration loop
if(verbose==TRUE){close(pb)}

  burn<-burn+1
	nCh                                     #number of chains
	Mu.burned<-Mus[[nCh]][burn:iterations]
	Cov.burned<-Covs[[nCh]][burn:iterations]
	Ps.burned<-Ps[[nCh]][burn:iterations]
	Zs.burned<-Zs[[nCh]][burn:iterations]
  Loglike.burned<-Loglike[burn:iterations]
	SteadyScore.burned<-SteadyScore$K0[burn:iterations]

# make weights and Zs matrices
Zs.burned<-t(sapply(Zs.burned, rbind))
Ps.burned<-t(sapply(Ps.burned, cbind))

	return(list(
    Mu = Mu.burned,
    Cov=Cov.burned,
    Ps= Ps.burned,
    Zs=Zs.burned,
    k.occupied=SteadyScore.burned,
    Log.likelihood=Loglike, y=y))
}
