
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
