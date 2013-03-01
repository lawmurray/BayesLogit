% package: BayesLogit
\name{compute.mixture}

\alias{compute.mixture}
\alias{compute.mixture.lg}
\alias{compute.mixture.l3}

\title{Compute Mixture}

\description{

  Compute the means, variances, and probabilities to approximate a
  -log(Ga) distribution or a type III logistic distribution.

}

\details{

  Fruhwirth-Schnatter et al. employ a table of discrete mixture of
  normals to approximate unruly families of distributions.  This
  function wraps to the code written by the Fruhwirth-Schnatter et al.
  to generate the means, variances, and probabilities that approximate a
  particular distribution.

  This function wraps to functions \code{compute.mixture.lg} and
  \code{compute.mixture.l3}.  \code{compute.mixture.l3} returns a list
  _without_ an array of means, since they are all identically zero.

}

\usage{

compute.mixture(shape, type=c("log.gamma", "logistic.iii"))

}

\arguments{

  \item{shape}{The shape parameter.}

  \item{type}{The type of distribution to approximate, either a
  -log(Ga(shape,1)) distribution or a type III logistic(shape)
  distribution.}
	
}

\value{

  Returns a list with components \code{m}, \code{v}, \code{p}
  representing the mean, variance, and probability of the specific
  discrete mixture of normals.

}

\examples{

## Approximate - log[Ga] using mixture of normals.
df = 9
nm = compute.mixture(df, "log.gamma")
nc = length(nm$m)

r  = sample.int(nc, 10, replace=TRUE, prob=nm$p)
e  = rnorm(10, nm$m[r], sqrt(nm$v[r]))

r.post = draw.indicators(e, nm)

## Approximate type III logistic using mixture of normals.
df = 9
nm = compute.mixture(df, "logistic.iii")
nc = length(nm$m)

r  = sample.int(nc, 10, replace=TRUE, prob=nm$p)
e  = rnorm(10, nm$m[r], sqrt(nm$v[r]))

r.post = draw.indicators(e, nm)

}

\seealso{
  \code{\link{draw.indicators}}
}

\references{

  Agnes Fussl.  binomlogit: Efficient MCMC for Binomial Logit Models.  (2012).

  Sylvia Fruhwirth-Schnatter. Finite Mixture and Markov Switch Models
  (2007), Springer.

  Sylvia Fruhwirth-Schnatter, Rudolf Fruhwirth, Leonhard Held, Havard
  Rue.  Improved auxiliary mixture sampling for hierarchical models of
  non-Gaussian data.  Statistics and Computing (2009).

  Sylvia Fruhwirth-Schnatter and Rudolf Fruhwirth. Data Augmentation
  and MCMC for Binary and Multinomial Logit Models.  Statistical
  Modelling and Regression Structures (2010), Springer-Verlag.

  Agnes Fussl, Sylvia Fruhwirth-Schnatter and Rudolf Fruhwirth.
  Efficient {MCMC} for Binomial Logit Models, (2013).

}

\keyword{normal mixture}