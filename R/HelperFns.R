#' GaussianLogPrior
#'
#' GaussianLogPrior is used to generate guassian prior of regression coefficients.
#'
#' @param x Regression Coefficients
#'
#' @param mu Mean
#'
#' @param sigma Standard Deviation
#'
#' @return GaussianLogPrior returns gaussian prior in log scale.
#'
#' @examples
#'  beta = runif(10)
#'  GaussianLogPrior(beta, 0, 1)
#'
#' @export
#'
#'
GaussianLogPrior <- function(x, mu, sigma){
  return(sum(-0.5*(x - mu)^2 / sigma^2))
}
