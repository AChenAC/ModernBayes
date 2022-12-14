% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Logistic_prediction.R
\name{Logistic.prediction}
\alias{Logistic.prediction}
\title{Logistic.prediction}
\usage{
Logistic.prediction(b, X, y)
}
\arguments{
\item{b}{beta samples returned by MCMC/SVGD.beta.generate.}

\item{X}{Selected X covariates returned by MCMC/SVGD.beta.generate.}

\item{y}{Outcome.}
}
\value{
MCMC.prediction returns predicted outcome with samples of beta coefficients.
}
\description{
Logistic.prediction is used to generate predictive outcomes of logistic regression.
}
