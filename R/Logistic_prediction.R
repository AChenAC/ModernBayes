#' Logistic.prediction
#'
#' Logistic.prediction is used to generate predictive outcomes of logistic regression.
#'
#' @param b beta samples returned by MCMC/SVGD.beta.generate.
#'
#' @param X Selected X covariates returned by MCMC/SVGD.beta.generate.
#'
#' @param y Outcome.
#'
#' @import stats
#'
#' @return MCMC.prediction returns predicted outcome with samples of beta coefficients.
#'
#' @export
#'
Logistic.prediction = function(b, X, y){
  pred = X %*% t(b)
  res_pred = apply(array(as.numeric((pred-y)>0), dim(pred)), 2, mean)
  return(factor(as.numeric(res_pred>0.5)))
}
