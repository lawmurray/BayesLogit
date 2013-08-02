## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

## This R script implements Fruhwirth-Schnatter and Fruhwirth's normal-mixture
## approximation to logistic regression (2007).  I mostly adhere to their
## notation.

## Define normal mixture.  FS&F p. 3511.
normal.mixture = list(
  w = c(0.00397, 0.0396, 0.168, 0.147, 0.125, 0.101, 0.104, 0.116, 0.107, 0.088),
  m = c(5.09, 3.29, 1.82, 1.24, 0.764, 0.391, 0.0431, -0.306, -0.673, -1.06),
  v = c(4.50, 2.02, 1.10, 0.422, 0.198, 0.107, 0.0778, 0.0766, 0.0947, 0.146)
  )
  normal.mixture$s = sqrt(normal.mixture$v)
  normal.mixture$N = length(normal.mixture$w)

ev.samp = -log(rexp(10000));
normal.mixture$mar.mean = mean(ev.samp);
normal.mixture$mar.var  = var(ev.samp);

## Make a copy.
NM = normal.mixture

draw.beta <- function(y.u, X, r, b.0=NULL, B.0=NULL, P.0=NULL)
{
  ## y: N x 1 outcomes.
  ## X: N x P design matrix.
  ## m: N x 1 means.
  ## s: N x 1 std. dev.
  ## b.0: prior mean for beta
  ## B.0: prior variance for beta
  ## P.0: prior precision for beta.
  
  ## FS-F use b to denote means and B to denote variances.

  N = nrow(X);
  P = ncol(X);

  if (is.null(b.0)) b.0 = rep(0.0, P);
  if (is.null(P.0)) P.0 = matrix(0.0, P, P);
  if (!is.null(B.0)) P.0 = solve(B.0);

  Xdv = X / NM$v[r];
  
  P.L = t(X) %*% Xdv;
  a.L = t(Xdv) %*% (y.u - NM$m[r]);

  P.N = P.0 + P.L;
  B.N = solve(P.N);
  b.N = B.N %*% (a.L + P.0 %*% b.0);
  
  beta = b.N + t(chol(B.N)) %*% rnorm(P)
} ## draw.beta

draw.utility <- function(y, lambda)
{
  ## y: N x 1 outcomes
  ## lambda = X beta
  N = length(y)
  U = runif(N); V = runif(N);
  y.u = -1 * log( -1 * log(U) / (1+lambda) - (1-y) * log(V) / lambda )
} ## draw.utility

draw.indicators.R <- function(y.u, lambda)
{
  ## y.u - N x 1 - latent variable y^u in paper.
  ## lambda = X beta

  res = y.u - log(lambda)
  log.wds = log(NM$w) - log(NM$s);

  ## Safer to work on log scale.  Columns correspond to outcome index i!
  log.post = -0.5 * outer(-1*NM$m, res, "+")^2 / NM$v + log.wds;
  unnrm.post = exp(log.post);

  ## Now sample.
  r = apply(unnrm.post, 2, function(x){sample.int(n=NM$N, size=1, prob=x)})
}  ## draw.indicators

logit.mix.gibbs <- function(y, X, samp=1000, burn=100, b.0=NULL, B.0=NULL, P.0=NULL, verbose=10000,
                            beta.true=NULL, y.u.true=NULL, r.true=NULL)
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
    y.u  = array(0, dim=c(M, N)),
    r    = array(0, dim=c(M, N))
    )

  beta = matrix(0.0, P);  # even odds.
  y.u  = matrix(0.0, N);
  r    = matrix(1.0, N);

  ## Always be careful with indicators.  You don't want to start with a very
  ## unlikely value (given the data).  It doesn't matter here because we can
  ## sample y.u without respect to r to start.  Nonetheless, seed randomly out
  ## of principle.
  r = sample.int(NM$N, N, prob=NM$w, replace=TRUE);

  ## In case we are doing testing.  May remove later.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(y.u.true))  y.u  = y.u.true;
  if (!is.null(r.true))    r    = r.true;
  
  for (i in 1:(samp+burn)) {

    lambda = exp(X %*% beta);

    ## (y.u, r | beta)  = (y.u | beta) (r | y.u, \beta)
    ## (y.u | beta, r) != (y.u | beta)
    y.u  = draw.utility(y, lambda)
    r    = draw.indicators.R(y.u, lambda)
    
    beta = draw.beta(y.u, X, r, b.0, "P.0"=P.0);
    
    if (i > burn) {
      out$beta[i-burn,] = beta;
      out$y.u [i-burn,] = y.u;
      out$r   [i-burn,] = r;
    }

    if (i %% verbose == 0) cat("Iteration", i, "\n");
    
  }

  out
}

################################################################################
                            ## Utility Functions ##
################################################################################

plot.compare <- function(out.fs, out.bl)
{
  P = ncol(out.fs$beta)
  old.par = par(mfrow=c(2,2));
  
  for (i in 1:P) {
    hist(out.fs$beta[,i], main="FS&F", xlab=paste("beta", i));
    acf(out.fs$beta[,i], main="FS&F");
    cat("mean:", mean(out.fs$beta[,i]), "sd:", sd(out.fs$beta[,i]), "\n");
    hist(out.bl$beta[,i], main="BL", xlab=paste("beta", i));
    acf(out.bl$beta[,i], main="BL");
    cat("mean:", mean(out.bl$beta[,i]), "sd:", sd(out.bl$beta[,i]), "\n");
    readline("Press <ENTER> to continue.");
  }

  par(mfrow=old.par);
}

## Original Plotting stuff.
## par(mfrow=c(2,2))
## hist(out.fs$beta[,1], main="FS&F")
## hist(out.fs$beta[,2], main="FS&F")
## hist(out.bl$beta[,1], main="BayesLogit")
## hist(out.bl$beta[,2], main="BayesLogit")

## readline("Press ENTER to continue.")

## par(mfrow=c(2,2))
## acf(out.fs$beta[,1], main="FS&F")
## acf(out.fs$beta[,2], main="FS&F")
## acf(out.bl$beta[,1], main="BayesLogit")
## acf(out.bl$beta[,2], main="BayesLogit")
             
################################################################################
                                   ## MAIN ##
################################################################################

## require("BayesLogit", lib.loc="../BLPackage/Test/")

## Synthetic test 1 ##

if (FALSE) {
  
  N = 1000;
  samp = 1000;
  burn = 500
  
  beta = c(1.0, 0.4);
  X = cbind(1, rnorm(N));

  r    = sample.int(NM$N, N, prob=NM$w, replace=TRUE);
  ep   = rnorm(N, NM$m[r], NM$s[r]);
  y.u  = X %*% beta + ep;
  y.u0 = -1 * log(rexp(N));
  y    = as.numeric(y.u > y.u0);
  
  B.0 = diag(1, 2)
  out.fs = logit.mix.gibbs(y, X, samp=samp, burn=burn, beta.true=beta, y.u.true=y.u, r.true=r);

  glm1   = glm(y ~ X+0, family=binomial(link="logit"));
  print(coef(glm1));
  out.bl = logit.gibbs(y, X, samp=samp, burn=burn, n.prior=0);

  plot.compare(out.fs, out.bl);
}

## Synthetic test 2 ##

if (FALSE) {
  
  N = 1000;
  
  beta = c(1.0, 0.4);
  X = cbind(1, rnorm(N));
  psi = X %*% beta;
  p = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, prob=p);

  out = logit.mix.gibbs(y, X, samp=1000, burn=0);

  glm1   = glm(y ~ X+0, family=binomial(link="logit"));
  ## out.bl = logit.gibbs(y, X, samp=1000, burn=0, n.prior=0);

}
