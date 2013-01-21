% package: BayesLogit
\name{CUBS}

\alias{CUBS.R}
\alias{CUBS.C}
\alias{FFBS.R}
\alias{FFBS.C}

\title{Conjugate Update Backward Sample}

\description{

  Draw an approximation of the states in a dyanmic generalized linear model
  using conjugate updating and a backwards sample.

}

\details{

  In a generalized linear model, the observation \eqn{y_i} is related,
  non-linearly, to some parameter \eqn{\psi_i}, which is modeled linearly:
  \eqn{x_i \beta}.  In a dynamic model, \eqn{beta} is allowed to change in time
  so that \eqn{\psi_i = x_i \beta_i} where \eqn{beta_i} is governed by some
  prespecified dynamics.  For the purposes of this routine we assume that
  \eqn{beta_i \sim AR(1)}.  One may approximately draw \eqn{\beta | y} by using
  linear Bayes with a conjugate updating procedure and the using backward
  sampling formulas.

}

\usage{

CUBS.C(z, X, n, mu, phi, W, m0, C0, obs=c("binom", "nbinom", "norm"), eps.rel=1e-8, max.iter=100)

}


\arguments{
  \item{z}{An N dimensional vector; \eqn{z_i} is the response at \eqn{x_i}.}
  \item{X}{An N x P dimensional design matrix.}
  \item{n}{An N dimensional vector; the interpreation of \eqn{n_i} changes with obs.}
  \item{mu}{A P.b dim. vector; the mean of beta.}
  \item{phi}{A P.b dim. vector; the autocorrelation of beta.}
  \item{W}{A P.b x P.b dim matrix; the innovation variance of beta.}
  \item{m0}{A P.b dimensional prior mean.}
  \item{C0}{A P.b x P.b dimensional prior precision.}
  \item{obs}{The type of response.}
  \item{eps.rel}{The relative error condition used when calculating the conjugate update parameters.}
  \item{max.iter}{The maximum number of iterations used when calculating the conjugate update parameters.}
}	

\value{

  \code{logit} returns a list.

  \item{beta}{A (N+1) x P array; the posterior sample of the dynamic regression
  coefficients.}

  \item{a.rate}{Acceptance rate.}

  \item{ess.time}{The time to generate the samples that are returned.}

}

\examples{

}

\references{
  Nicholas G. Polson, James G. Scott, and Jesse Windle.  Bayesian inference for
  logistic models using Polya-Gamma latent variables.
  \url{http://arxiv.org/abs/1205.0310}

  Nicholas Poslon and James G. Scott.  Default Bayesian analysis for
  multi-way tables: a data-augmentation approach.
  \url{http://arxiv.org/pdf/1109.4180}
}

\keyword{logit}
\keyword{regression}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.
