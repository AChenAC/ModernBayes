#' SVGD.beta.generate
#'
#' SVGD.beta.generate is used to generate samples of beta coefficients for top 1/2 features selected using PCA.
#'
#' @param model Regression model.
#'
#' @param prior Prior functions.
#'
#' @param likelihood Log likelihood function.
#'
#' @param X Sets of predictors.
#'
#' @param y Outcome.
#'
#' @import randtoolbox
#'
#' @return MCMC.beta.generate returns the samples of beta coefficients.
#'
#' @export
#'
SVGD.beta.generate = function(model, prior, likelihood, X, y){
  pca = prcomp(X, center = TRUE,scale. = TRUE)
  X = pca$x
  X = cbind(1, X)[,1:(ncol(X)/2+1)]
  d = dim(X)[2]
  n = length(y)
  n_particles = 30

  sigmoid = function(s) {
    return(1 / (1 + exp(-s)))
  }

  dsigmoid = function(s) {
    return(1 - sigmoid(s))
  }

  prior_mu = 0.
  prior_sigma = 1.

  gradLogPosteriorFunction = function(beta, X, y){ # d/db log Pr(b|X,y), specific to logistic likelihood and normal prior
    y_samples = array(rep(y, n_particles), c(n, n_particles))
    z = X %*% t(beta)
    y_pred_samples = sigmoid(z)
    grad_outer = (y_samples-y_pred_samples) / (y_pred_samples-y_pred_samples^2 + 1e-3)
    grad_posterior = t(dsigmoid(z) * grad_outer) %*% X / n + (prior_mu - beta) / prior_sigma^2
    return(grad_posterior)
  }

  get_gradient = function(beta, X, y){ # svgd gradient
    grad_posterior = gradLogPosteriorFunction(beta, X, y)
    Q = t(grad_posterior) %*% grad_posterior / n_particles

    sign_diff = array(NA, c(n_particles, n_particles, d))
    for(di in 1:d){
      sign_diff[,,di] = outer(beta[,di], beta[,di], '-')
    }
    Qdiff = array(array(sign_diff, c(n_particles^2, d)) %*% Q, c(n_particles, n_particles, d))
    diffQdiff = sign_diff*Qdiff
    diffQdiffsum = array(NA, c(n_particles, n_particles))
    for(ni in 1:n_particles){
      diffQdiffsum[ni,] = apply(diffQdiff[ni,,], 1, sum)
    }
    h = mean(diffQdiffsum) / log(n_particles)

    kern  = exp(-diffQdiffsum / (2*h))
    gkern = array(NA,c(n_particles, d))
    for(ni in 1:n_particles){
      gkern[ni,] = -t(Qdiff[,ni,]) %*% kern[ni,]/ h
    }
    mgJ = kern %*% grad_posterior + gkern
    return(mgJ / n_particles)
  }

  alpha = 1. # gradient descent step size
  tol = 1e-5

  beta_samples = prior_mu + sobol(n_particles, d, normal=TRUE)*prior_sigma # initial beta particles sampled from prior
  adag = array(0, dim=dim(beta_samples)) # adagrad parameter

  db = c()
  for(i in 1:10000){
    delta_beta = get_gradient(beta_samples, X, y) # svgd gradient
    adag = adag + delta_beta^2 # adagrad parameter
    beta_samples = beta_samples + alpha*delta_beta / sqrt(adag) #gradient descent via adagrad
    db_i = norm(delta_beta) / (n_particles*d) # norm of gradient, for convergence criteria
    db = c(db, db_i)
    if(db_i < tol){ # convergence criteria
      break
    }
  }
  return(list(X=X, beta_samples = beta_samples))
}
