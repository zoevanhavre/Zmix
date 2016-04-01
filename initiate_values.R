ZZ<-init.Z


initiate_values<-function(ZZ){
  IndiZ <-  (ZZ == matrix((1:k), nrow = n, ncol = k, byrow = T))
  ns <-     apply(IndiZ,2,sum)  	#  size of each group

  .Ysplit<- replicate(k, list())	#storage to group Y's
  WkZ<-	    replicate(k, list(NA))	#storage for within group variability
  ybar<-	  replicate(k, list(NA))

  if(r==1){
  sx <- apply(IndiZ*y, 2, sum) } else{ sx<-NA}

  for (r>1 && .i in 1:k){
    .Ysplit[[.i]]<-y[ZZ==.i,]
    if (ns[.i]>1){					# for groups with >1 obsevations
      ybar[[.i]]<- as.matrix(t(apply(.Ysplit[[.i]], 2, mean)))
      } else if (ns[.i]==1){
      ybar[[.i]]<-  t( as.matrix(.Ysplit[[.i]]))
    }
    #Within group unexplained variability
    if (ns[.i]==1){
      WkZ[[.i]]<-crossprod(as.matrix(.Ysplit[[.i]]-ybar[[.i]]))
      } else if (ns[.i]>1){
      WkZ[[.i]]<-0
        for (.n in 1:ns[.i]){
        WkZ[[.i]]<-WkZ[[.i]]+ crossprod( .Ysplit[[.i]][.n,]-ybar[[.i]])
        }
      } else if (ns[.i]==0){
      WkZ[[.i]]<-NA
    }

  }
      list('ns'=ns,'ybar'=ybar, 'WkZ'=WkZ, 'sx'=sx)
  }
