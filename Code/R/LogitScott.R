## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.




normal.llh <- function(x, m, S)
{
  ## n = length(x)
  ## evd = eigen(S, symmetric=TRUE);
  ## y = t(evd$vectors) %*% (x - m);
  ## llh = -0.5 * n * log(2*pi) - 0.5 * sum(log(evd$values)) - 0.5 * t(y) %*% (y / evd$values);
  llh = dmvnorm(as.numeric(x), as.numeric(m), S, log=TRUE)
  llh
}

ind.logistic.llh <- function(x, m, s)
{
  ## n = length(x);
  ## z = (x - m) / s;
  ## llh = -2.0 * sum( log(cosh(0.5 * z)) ) - 2 * sum(log(s)) + 2 * n * log(0.5);
  llh = sum(dlogis(x, m, s, log=TRUE));
  llh
}

draw.z <- function(lambda, y){
  n = length(lambda)
  u = runif(n)
  z = log(lambda * u + y) - log(1 - u + lambda * (1-y));
  z
} ## draw.z

logit.DFAE.gibbs <- function(y, X, b.0=NULL, B.0=NULL, P.0=NULL, samp=1000, burn=500, verbose=500,
                             beta.true=NULL, z.true=NULL)
{
  N = nrow(X)
  P = ncol(X)
  M = samp

  ## Default prior parameters.
  if (is.null(b.0)) b.0 = rep(0.0, P);
  if (is.null(P.0)) P.0 = matrix(0.0, P, P);
  if (!is.null(B.0)) P.0 = solve(B.0);
  
  out = list(
    beta = array(0, dim=c(M, P)),
    z    = array(0, dim=c(M, N)),
    accept = rep(0, M)
    )

  beta = matrix(0.0, P);  # even odds.
  z    = matrix(0.0, N);

  num.accept = 0;
  accept = 0;

  ## In case we are doing testing.  May remove later.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(z.true))    z = z.true;
  
  ## Preprocess
  XX = t(X) %*% X
  sig2 = pi^2 / 3;
  ## C.0 = solve(P.0);
  P.L  = XX / sig2;
  P.N = P.L + P.0;
  evd = eigen(P.N);
  S = evd$vectors %*% diag(1.0 / evd$values) %*% t(evd$vectors);
  ## S = solve(P.N);
  U   = chol(S)

  for (i in 1:(burn+samp)) {

    psi = X %*% beta;
    lambda = exp(psi);

    ## Draw z.
    z = draw.z(lambda, y);

    ## Calculate LS+Prior estimate.
    m      = S %*% (t(X) %*% z / sig2 + P.0 %*% b.0)
    beta.p = m + t(U) %*% rnorm(P);
    psi.p  = 1 * X %*% beta.p

    ## Calculate acceptance prob.
    ## log.prop.ratio = ind.logistic.llh(z, psi.p, 1.0)+normal.llh(beta.p, b.0, C.0)-normal.llh(beta.p, m, S);
    ## log.prev.ratio = ind.logistic.llh(z, psi  , 1.0)+normal.llh(beta  , b.0, C.0)-normal.llh(beta  , m, S);

    diff1 = ind.logistic.llh(z, psi.p, 1.0) - ind.logistic.llh(z, psi  , 1.0)
    diff2 = -0.5 * t(beta.p - b.0) %*% P.0 %*% (beta.p - b.0) + 0.5 * t(beta - b.0) %*% P.0 %*% (beta - b.0)
    diff3 = -0.5 * t(beta - m) %*% P.N %*% (beta - m) + 0.5 * t(beta.p - m) %*% P.N %*% (beta.p - m)

    ratio = exp(diff1 + diff2 + diff3);
    alpha = min(1, ratio);

    accept = 0
    if (runif(1) < alpha) {
      beta = beta.p
      num.accept = num.accept + 1
      accept = 1
    }

    if (i > burn) {
      out$beta[i-burn,] = beta;
      out$z[i-burn, ]  = z;
      out$accept[i-burn] = accept;
    }

    if (i %% verbose == 0) cat("Iteration:", i, "Accept Rate:", num.accept/i, ".\n");
    
  }

  out$rate = mean(out$accept);
  
  out
}

################################################################################

if (FALSE) {
  
  N = 1000;
  samp = 1000;
  burn = 500;

  P.0 = diag(0.01, 2);
  
  beta = c(1.0, 0.4);
  X = cbind(1, rnorm(N));
  psi = rks(N);
  lambda = (2*psi)^2;
  ep = rnorm(N, 0.0, sqrt(lambda));
  z  = X %*% beta + ep
  y  = as.numeric(z > 0);

  start = proc.time()
  out = logit.DFAE.gibbs(y, X, samp=samp, burn=burn, verbose=500, P.0=P.0)
  total.time = proc.time() - start;

  glm1   = glm(y ~ X+0, family=binomial(link="logit"));
  summary(glm1)
  ## out.bl = logit.gibbs(y, X, samp=samp, burn=burn, n.prior=0);

}

if (FALSE) {
  
  N = 500;
  samp = 5000;
  burn = 2000;
  
  beta = c(1.0, 0.4, 1.2, 0.8, -1.0, 0.5, -1.5);
  X = cbind(1, matrix(rnorm(N*6), ncol=6));
  psi = X %*% beta;
  ep = rlogis(N, 0, 1)
  z  = X %*% beta + ep
  y  = as.numeric(z > 0);

  P.0 = diag(0.01, 7);
  ## P.0 = diag(0.00, 7);
  
  start = proc.time()
  out = logit.DFAE.gibbs(y, X, samp=samp, burn=burn, verbose=500, P.0=P.0)
  ## out = logit.DFAE.gibbs(y, X, samp=samp, burn=burn, verbose=500, P.0=P.0, z.known=z.mean)
  total.time = proc.time() - start;

  glm1   = glm(y ~ X+0, family=binomial(link="logit"));
  summary(glm1)
  ## out.bl = logit.gibbs(y, X, samp=samp, burn=burn, n.prior=0);

}

################################################################################



################################################################################
                                 ## APPENDIX ##
################################################################################

## Initially, I had misspecified the densities.  That throws everything way off.
