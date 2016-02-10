#' Processing output 
#'
#' inner function.
#' @param x: vector of factor values  (here, allocations)
#' @keywords postprocessing
#' @export
#' @examples
#' #NA

			maxZ<-function (x) {
			as.numeric(names(which.max(table( x ))))	
			} 
