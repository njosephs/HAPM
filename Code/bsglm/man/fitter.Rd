% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsglm.R
\name{fitter}
\alias{fitter}
\title{Fit a Bayesian generalized linear model}
\usage{
fitter(
  dispersion = 1,
  adjust_dispersion = FALSE,
  prior_coef = NULL,
  prior_disp = NULL
)
}
\arguments{
\item{dispersion}{value of dispersion if assumed deterministic or initial
value for fitting}

\item{adjust_dispersion}{is dispersion random?}

\item{prior_coef}{hyper-prior parameters for the coefficients, specified as
a list with entries \code{mean}, a vector, and \code{precision}, either a matrix or a
vector storing the diagonal entries}

\item{prior_disp}{hyper-prior parameters for the dispersion, specified as a
list with entries \code{df}, the number of degrees of freedom, and \code{scale}, the
variance prior scale}
}
\value{
a \code{glm} object.
}
\description{
Find the maximum a posteriori (MAP) estimator of Bayesian generalized
linear model with coefficients \eqn{\beta}{beta} and dispersion
\eqn{\phi}{\phi}. The coefficients have a normal prior with mean
\eqn{\beta_0}{beta0} and precision matrix \eqn{\Omega}{Omega},
\eqn{\beta \sim N(\beta_0, \Omega^{-1})}{beta ~ N(beta0, Omega^{-1})}, and
the dispersion has an inverse scaled chi-square prior with degrees of
freedom \eqn{\nu}{nu} and scale \eqn{\tau^2}{tau2},
\eqn{1/\phi \sim Ga(\nu/2, \nu\tau^2/2)}{1/phi ~ Ga(nu / 2, nu * tau2 / 2)}.
}
\details{
This routine is a \code{glm} method and should be passed to a \code{glm} call, e.g.,
\code{glm(formula, data, family, ..., method = bsglm::fitter(...))}.
}
