#' MCMC.beta.generate
#'
#' MCMC.beta.generate is used to generate samples of beta coefficients for top half features selected using PCA.
#'
#' @param model Regression model.
#'
#' @param prior Prior functions.
#'
#' @param likelihood Log likelihood function.
#'
#' @param posterior Posterior function.
#'
#' @param X Sets of predictors.
#'
#' @param y Outcome.
#'
#' @return MCMC.beta.generate returns the samples of beta coefficients.
#'
#' @export
#'
MCMC.beta.generate = function(model, prior, likelihood, posterior, X, y){
  # #data preprocessing
  X = scale(X) # centered features
  pca = prcomp(X, center = TRUE,scale. = TRUE) # pca on centered features (not essential, but helpful)
  X = pca$x[,1:(ncol(X)/2)] # top half principle component features
  d = dim(X)[2]
  #sample initialization
  n = 100000 # number of mcmc samples
  n_burnin = n / 10 # 10% burn-in
  samples = matrix(0, nrow=(n+n_burnin), ncol=d+1) # allocation
  currentProbability = posterior(samples[1,], X, y) # mcmc p(b | X, y) sample 1
  proposal_cov_chol = chol(diag(d+1))
  sd = 2.4^2 / d # scaling factor, see Haario et al 2001

  for(i in 1:(n + n_burnin-1)){ # vanilla MCMC

    if(i==n_burnin){
      proposal_cov_chol = chol(cov(samples[1:n_burnin,])*sd)
    }
    deltaProposal = proposal_cov_chol %*% rnorm(d+1)

    proposal = matrix(samples[i,] + deltaProposal, nrow=d+1)
    proposalProbability = posterior(proposal, X, y)
    logAlpha = (proposalProbability - currentProbability) # log alpha = log proposalP - log currentP

    logU = log(runif(1)) # log u
    acceptance = logAlpha >= logU # bool

    if(is.na(acceptance)){ # if nan
      samples[i+1,] = rnorm(d+1, 0, 1) # if nan then start chain at random value
      currentProbability = posterior(samples[i+1,], X, y)
    }

    else if(acceptance){ # if accept
      # print(TRUE)
      samples[i+1,] = proposal
      currentProbability = proposalProbability
    }

    else if(!acceptance){ # if reject
      # print(FALSE)
      samples[i+1,] = samples[i,]
    }
  }
  sample_stationary = samples[n_burnin:(n+n_burnin),]
  return(list(X = cbind(1, X), sample_stationary=sample_stationary))
}
