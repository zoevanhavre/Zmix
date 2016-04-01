
# univariate

## univ priors
#tau=     1
n<-       length(y)
a<-       2.5
b<-       2/var(y)
d<-       2
lambda<-  sum(y)/n

# initiate
mux<-list(mu=seq(from=min(y), to=max(y),length.out=k),sigma=rep(1, k),p=rep(1/k,k), k=k)
mu0<- mux$mu
p0<-  mux$p
sig0<-mux$sigma

  for (.ch in 1:nCh){
  Bigmu[[.ch]][1,] <-     mux$mu                # initial value of mu's
  Bigp[[.ch]][1,] <-      mux$p                 # inital value of p's
  Bigsigma[[.ch]][1,] <-  mux$sigma             # inits for sigma
	}

#needed: ns, sx, sv
sx <-   apply(IndiZ*y, 2, sum)
.bmu<-  matrix((1:k), nrow = n, ncol = k, byrow = T)
for (t in 1:n) {
        .bmu[t,]<-Bigmu[[.ch]][j,]  }
sv<-apply((y*IndiZ-.bmu*IndiZ)^2, 2, sum)

### main FUNCTIONS

# pzs
# Iter=1
for (i in 1:n) {
Pzs[[.ch]][i,] <- (p0/sqrt(sig0))*exp(-((Y[i]-mu0)^2)/(2*sig0))
Pzs[[.ch]][i,] <- Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,]) }											  # Scale to equal 1?
# iter>1
for (i in 1:n) {
Pzs[[.ch]][i,] <- (Bigp[[.ch]][j-1,]/sqrt(Bigsigma[[.ch]][j-1,]))*exp(-((Y[i]-Bigmu[[.ch]][j-1,])^2)/(2*Bigsigma[[.ch]][j-1,]))                                                                                                    # Scale to equal 1
Pzs[[.ch]][i,] <- Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,])
}


# other pars
Bigp[[.ch]][j,]<-     rdirichlet(m=1,par= ns+alphas[.ch])
Bigp[[.ch]][j,] <-    rdirichlet(m=1,par= ns+alphas[.ch])
# Bigsigma[[.ch]][j-1,] , or sig0 if iter=1
Bigmu[[.ch]][j,]<-    rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][j-1,]/(tau+ns)))
Bigsigma[[.ch]][j,]<- rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)

# update for iter>1
for (i in 1:length(Bigmu[[.ch]][j,])){ if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }
.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)                                                         # changes, added /ns

Bigsigma[[.ch]][j,]<- rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)
