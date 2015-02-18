#' Draw from Inverse Gamma 
#'
#' This function draws samples from an inverse Gamma dist
#' @param m Vector of Alpha values (hyperparameter of dist)
#' @keywords Gamma
#' @export
#' @examples
#' rinvgamma(1, 2)


 rinvgamma<-function (n, shape, scale = 1) {
		  		  return(1/rgamma(n = n, shape = shape, rate = scale))}

