	 
#' All permutations given size 
#'
#' This function outputs all permutations of a set number of values
#' @param n, number of values in permutation
#' @keywords permutationa
#' @export
#' @examples
#' permutations(2)

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