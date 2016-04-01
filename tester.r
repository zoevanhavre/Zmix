


source("/Users/vanhavre/Documents/Zmix/funks.R")
source("/Users/vanhavre/Documents/ZmixFix.R")

#####
y<-iris[,1:2]



iterations=100
k=10
alphas= c(1/2^(c(6, 10, 30)))
burn=50
init.method="Kmeans"
verbose=TRUE



# test the thing


wrapFunk<-function(x, k=10){
  innerFunk(x)
  }

  innerFunk<-function(a){
    a+k
  }

wrapFunk(1)

#Error in innerFunk(x) : object 'k' not found

innerFunk<-function(a, k){a+k}

wrapFunk<-function(x, k=10){
  innerFunk(x, k)
  }

wrapFunk(1)



# using dot dot dot
tf<-function(...){
list(...) ->pars
  return(pars )
}

tf(1,2,3)



innerFunk<-function(...){
a<-list(...)[[1]]
k<-list(...)[[2]]

a+k}

wrapFunk<-function(x, k=10){
  innerFunk(x, k)
  }

wrapFunk(1)

# OK, works
