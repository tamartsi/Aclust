#' Title expit
#' 
#' Returns inverse of the logistic link function
#'
#' @param x numeric vector
#'
#' @return exp(x)/(1+exp(x))
#' @export
#'
#' @examples
expit <-
function (x) 
{
    exp(x)/(1 + exp(x))
}
