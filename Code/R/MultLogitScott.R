


norm.llh <- function(x, m, S)
{
  n = length(x)
  evd = eigen(S, symmetric=TRUE);
  y = t(evd$vectors) %*% x;
  llh = -n * log(2*pi) - 0.5 * sum(log(evd$values)) - 0.5 * t(y) %*% (y / evd$values);
  llh
}

iid.logistic.llh <- function(x, m, s)
{
  n = length(x);
  z = (x - m) / s;
  llh = -1.0 * sum( log(cosh(0.5 * z)) ) - n * log(s);
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
    z    = array(0, dim=c(M, N))
    )

  beta = matrix(0.0, P);  # even odds.
  z    = matrix(0.0, N);

  ## In case we are doing testing.  May remove later.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(z.true))    z = z.true;
  
  ## Preprocess
  XX = t(X) %*% X
  sig2 = pi^2 / 3;
  P.L  = XX / sig2;
  reg.mat = solve(XX) %*% t(X);
  P.N = P.L + P.0;
  S   = solve(P.N);
  C   = chol(S)
  
  for (i in 1:(burn+samp)) {

    psi = X %*% beta;
    lambda = exp(psi);

    ## Draw z.
    z = draw.z(lambda, y);

    ## Calculate LS+Prior estimate.
    m = S %*% (t(X) %*% z / sig2 + P.0 %*% b.0)
    beta.p = m + C %*% rnorm(P);
    psi.p = X %*% beta.p

    ## Calculate acceptance prob.
    log.logis.ratio = iid.logistic.llh(z, psi.p, 1.0) - iid.logistic.llh(z, psi, 1.0);
    log.prior.ratio = 0.0
    log.q.ratio     = normal.llh(z, psi.p, S) - normal.llh(z, psi, S);

    ratio = exp(log.logis.ratio - log.q.ratio);
    alpha = min(1, ratio);

    if (runif(1) < alpha) beta = beta.p

    if (i > burn) {
      out$beta[i-burn,] = beta;
      out$z[i,-burn, ]  = z;
    }

    if (i %% verbose == 0) cat("Iteration", i, ".\n");
    
  }
  
}
