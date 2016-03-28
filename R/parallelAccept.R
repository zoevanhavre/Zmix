#' Prior Parallel Tempering: compute acceptance probabilities
#'
#' This function takes the mixture weights of two chains at one iteration, and the values of the hyperparameter alpha in these chains, and outputs the acceptance probability required for prior parallel tempering.
#' @param w1, w2, a1, a2
#' @keywords prior parallel tempering.
#' @export
#' @examples
#' # Coming zoon

parallelAccept<-function(w1, w2, a1, a2){
    w1[w1< 1e-200]<-1e-200             # truncate very small weights to be non-zero
    w2[w2< 1e-200]<-1e-200
    T1<-dDirichlet(w2, a1, log=TRUE)
    T2<-dDirichlet(w1, a2, log=TRUE)
    B1<-dDirichlet(w1, a1, log=TRUE)
    B2<-dDirichlet(w2, a2, log=TRUE)
    MH<-min(1,	exp(T1+T2-B1-B2))
    Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
    return(Ax)}
