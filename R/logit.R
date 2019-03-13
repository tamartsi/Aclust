#' Title logit
#'
#' Compute generalized logit function
#' 
#' @param x value(s) to be transformed 
#'
#' @return
#' @export
#'
#' @examples
logit <-
function (x) 
{
    log(x/(1 - x))
}
