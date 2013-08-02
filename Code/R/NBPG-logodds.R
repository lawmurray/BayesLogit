## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


################################################################################
                            ## When Modeling PSI ##
################################################################################

df.llh <- function(d, psi, G, ymax)
{
  p =  1 / (1 + exp(-psi))                                
  sum( log(d+0:(ymax-1)) * G[1:ymax] ) + d * sum(log(1-p));
}

draw.df <- function(d.prev, psi, G, ymax)
{
  d.new = d.prev

  ## Kernel 1: RW MH.
  ## And a random walk.
  rw.lower = max(d.prev - 1, 1);
  rw.upper = d.prev + 1;
  rw.grid  = rw.lower:rw.upper;
  rw.n     = length(rw.grid)
  rw.p     = rep(1/rw.n, rw.n);
  
  d.prop = sample(rw.grid, 1, prob=rw.p);
  
  ltarget = df.llh(d.prop, psi, G, ymax) - df.llh(d.prev, psi, G, ymax)
  lpps = log(ifelse(d.prop==1, 1/2, 1/3)) - log(ifelse(d.prev==1, 1/2, 1/3));
  lratio = ltarget + lpps

  if (runif(1) < exp(lratio)) d.new = d.prop
  
  d.new
}

draw.df.real <- function(Y, r, Pr) {
  ## r = d
  if(r > 1) rstar = runif(1,r-1,r+1)
  else rstar = runif(1,0,2)
  ll = sum(dnbinom(Y, r, Pr, log=TRUE))
  llstar = sum(dnbinom(Y, rstar, Pr, log=TRUE))
  lalpha = llstar - ll
  lu = log(runif(1))
  if(lu < lalpha) r = rstar
  r
}

################################################################################

draw.beta <- function(X, kappa, w, b.0=NULL, B.0=NULL, P.0=NULL)
{
  ## X: N x P design matrix.
  ## b.0: prior mean for beta
  ## B.0: prior variance for beta
  ## P.0: prior precision for beta.
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(b.0)) b.0 = rep(0.0, P);
  if (is.null(P.0)) P.0 = matrix(0.0, P, P);
  if (!is.null(B.0)) P.0 = solve(B.0);

  P.N = t(X) %*% (X * w) + P.0;
  a.L = t(X) %*% kappa;
  ## S = solve(PC); ## chol2inv works better for larger P?
  S = chol2inv(chol(P.N));
  m = S %*% (a.L + P.0 %*% b.0);
  beta = m + t(chol(S)) %*% rnorm(P);

} ## draw.beta

NB.PG.gibbs <- function(y, X,
                        b.0=NULL, B.0 = NULL, P.0 = NULL,
                        samp=1000, burn=500, verbose=500,
                        beta.true = NULL, w.true = NULL, d.true=NULL)
{
  ## X: n by p matrix
  ## y: n by 1 vector, counts.

  X = as.matrix(X);
  y = as.matrix(y);
  
  P = ncol(X)
  N = nrow(X)
  
  ## Default prior parameters.
  if (is.null(b.0)) b.0 = rep(0.0, P);
  if (is.null(P.0)) P.0 = matrix(0.0, P, P);
  if (!is.null(B.0)) P.0 = solve(B.0);
  
  ## Initialize.
  w = rep(0,N)
  beta = rep(0.0, P)
  d = 1;

  ## Set known.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(w.true))    w    = w.true;
  if (!is.null(d.true))    d    = d.true;

  ## Preprocess
  ymax = max(y);
  F = cumsum(hist(y, breaks=0:(ymax+1)-0.5, plot=FALSE)$counts)
  G = N - F;
  
  output <- list(w = matrix(nrow=samp, ncol=N),
                 beta = matrix(nrow=samp, ncol=P),
                 d = rep(0, samp)
                 )

  ## Sample
  for ( j in 1:(samp+burn) )
  {
    ## WARNING: JOINT DRAW.
    ## draw (d, w | beta)
    ## draw (d | beta)
    psi = drop(X%*%beta)
    pr  = 1 / (1 + exp(-psi))
    d = draw.df(d, psi, G, ymax);
    ## d = draw.df.real(y, d, pr)
    ## draw (w | d, beta)
    w = rpg.devroye(N, y+d, psi);

    ## draw beta
    kappa = 0.5 * (y-d)
    beta = draw.beta(X, kappa, w, b.0=b.0, P.0=P.0);
    
    # Record if we are past burn-in.
    if (j>burn) {
        output$w[j-burn,]    = w
        output$beta[j-burn,] = beta
        output$d[j-burn]     = d
    }

    if (j %% verbose == 0) { print(paste("NBPG-logodds: Iteration", j)); }
  }

  output
} ## logit.PG.gibbs

################################################################################

if (FALSE) {

  library("BayesLogit")
  ## source("NBPG-logodds.R")
  
  N = 500
  X = cbind(1, rnorm(N))

  beta = c(2.0, 1.5);

  ## Poisson-Gamma mixture.
  psi = X %*% beta;
  d  = 4
  alpha = exp(psi)
  lambda = rgamma(N, d, scale=alpha);
  y = rpois(N, lambda);

  samp = 1000
  burn = 500
  verbose = 500
  
  out <- NB.PG.gibbs(y, X,
                     samp=samp, burn=burn, verbose=verbose,
                     beta.true=NULL, w.true=NULL, d.true=NULL)
  
}


################################################################################
                                 ## APPENDIX ##
################################################################################

  ## if (k==3) {
  ##   ## Kernel 3: This is not a MH kernel since it won't traverse the
  ##   ## space.  This is kind of too complicated -- a discrete normal
  ##   ## approximation.
  ##   d.m = max(mle, 1)
  ##   d.v = fim
  ##   d.s = max(sqrt(d.v), 1/3);
  ##   d.lower = max(floor(d.m - 3 * d.s), 1);
  ##   d.upper = ceil(d.m + 3 * d.s);
  ##   d.grid  = d.lower:d.upper;
  ##   d.p     = dnorm(d.grid, d.m, d.s);
  ##   d.p     = d.p / sum(d.p);

  ##   d.prop = sample(d.grid, 1, prob=d.p);

  ##   if (d.prev < d.lower || d.prev > d.upper) {
  ##     d.new = d.prop
  ##   } else {
  ##     p.prop = 1 / (1 + d.prop / mu)
  ##     p.prev = 1 / (1 + d.prev / mu)
  ##     ltarget = dnbinom(y, d.prop, 1-p.prop, log=TRUE) - dnbinom(y, d.prev, 1-p.prev, log=TRUE)
  ##     lppsl   = dpois(d.prev, d.m, log=TRUE) - dpois(d.prop, d.m, log=TRUE)
  ##     lratio  = ltarget + lppsl
  ##   }
  ## }
