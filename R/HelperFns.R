#' GaussianLogPrior
#'
#' GaussianLogPrior is used to generate guassian prior of regression coefficients.
#'
#' @param x Regression Coefficients
#'
#' @param mu Mean. The default value is 0 if not specified.
#'
#' @param sigma Standard Deviation. The default value is 1 if not specified.
#'
#' @return GaussianLogPrior returns gaussian prior in log scale.
#'
#' @examples
#'  beta = runif(10)
#'  GaussianLogPrior(beta, 0, 1)
#'
#' @export
#'
GaussianLogPrior <- function(x, mu = 0, sigma = 1){
  return(sum(-0.5*(x - mu)^2 / sigma^2))
}


#' LaplaceLogPrior
#'
#' LaplaceLogPrior is used to generate Laplace prior of regression coefficients.
#'
#' @param x Regression Coefficients
#'
#' @param mu location parameter
#'
#' @param b scale parameter
#'
#' @return LaplaceLogPrior returns Laplace prior in log scale.
#'
#' @examples
#'  beta = runif(10)
#'  LaplaceLogPrior(beta, 0, 1)
#'
#' @export
#'
LaplaceLogPrior <- function(x, mu, b){
  return(sum(-abs(x - mu)/b))
}

################################################################################

#' LogisticModel
#'
#' LogisticModel is the mean function (also known as sigmoid function) of logistic regression.
#'
#' @param X Sets of predictors in matrix form.
#'
#' @param b Beta Coefficients in matrix form.
#'
#' @return LogisticModel returns the log version of the mean function of logistic regression.
#'
#' @export
#'
LogisticModel <- function(b, X){
  et = exp(-(cbind(1, X) %*% b))
  return(log(1 / (1 + et)))
}

#' GaussianModel
#'
#' GaussianModel is the mean function of gaussian regression.
#'
#' @param X Sets of predictors in matrix form.
#'
#' @param b Beta Coefficients in matrix form.
#'
#' @return GaussianModel returns the mean function of gaussian regression.
#'
#' @export
#'
GaussianModel <- function(b, X){
  return(cbind(1, X) %*% b)
}

#' GammaModel
#'
#' GammaModel is the mean function of gamma regression.
#'
#' @param X Sets of predictors in matrix form.
#'
#' @param b Beta Coefficients in matrix form.
#'
#' @return GammaModel returns the mean function of gamma regression.
#'
#' @export
#'
GammaModel <- function(b, X){
  return(1/(cbind(1, X) %*% b))
}

################################################################################
#' LogisticLogLik
#'
#' LogisticLogLik produces logistic probability.
#'
#' @param X Sets of predictors in matrix form.
#'
#' @param y outcome variable in matrix form.
#'
#' @param b Beta Coefficients in matrix form.
#'
#' @return LogisticLogLik returns the logistic log likelihood.
#'
#' @export
#'
LogisticLogLik <- function(b, X, y) { # logistic probability p(X, y | b)
  return((sum(LogisticModel(b, X[y,])) + sum(1-LogisticModel(b, X[!y,]))))
}

#' GaussianLogLik
#'
#' GaussianLogLik produces gaussian probability.
#'
#' @param X Sets of predictors in matrix form.
#'
#' @param y outcome variable in matrix form.
#'
#' @param b Beta Coefficients in matrix form.
#'
#' @param sigma standard deviation. The default value is 1 if not specified.
#'
#' @return LogisticLogLik returns the logistic log likelihood.
#'
#' @export
#'
GaussianLogLik <- function(b, X, y, sigma = 1){
  n = length(y)
  return(-n/2*log(2*pi) - n/2*log(sigma) - 1/2/sigma*sum((y - GaussianModel(b, X))^2))
}

################################################################################

#' logPosteriorFunction
#'
#' logPosteriorFunction produces un-normalized posterior pdf based on gaussian prior.
#'
#' @param X Sets of predictors in matrix form.
#'
#' @param y outcome variable in matrix form.
#'
#' @param b Beta Coefficients in matrix form.
#'
#' @return logPosteriorFunction returns un-normalized posterior pdf.
#'
#' @export
#'
logPosteriorFunction <- function(b, X, y){ # un-normalized posterior pdf
  return(logPriorFunction(b, X, y) + logLikelihoodFunction(b, X, y)) # return scalar p(b | X, y)
}
