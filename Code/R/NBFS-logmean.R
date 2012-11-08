
## Negative Binomial Regression using discrete mixture of normals.

## if (exists("TESTING")) {
source("ComputeMixture.R")
source("NB-Indicators.R");
source("NB-Shape.R") ## Routine for sampling shape.
## } ## TESTING

################################################################################
                             ## When Modeling MU ##
################################################################################

##------------------------------------------------------------------------------

draw.beta <- function(X, lambda, d, r, nmix, b.0=NULL, B.0=NULL, P.0=NULL)
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

  z = log(lambda) + log(d) + nmix$m[r];
  P.N = t(X) %*% (X / nmix$v[r]) + P.0;
  a.L = t(X) %*% (z / nmix$v[r]);
  ## S = solve(PC); ## chol2inv works better for larger P?
  S = chol2inv(chol(P.N));
  m = S %*% (a.L + P.0 %*% b.0);
  beta = m + t(chol(S)) %*% rnorm(P);

} ## draw.beta

NB.FS.gibbs <- function(y, X,
                        b.0=NULL, B.0 = NULL, P.0 = NULL,
                        samp=1000, burn=500, verbose=500,
                        beta.true = NULL, d.true=NULL, lambda.true=NULL, r.true=NULL)
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
  beta = rep(0.0, P)
  d = 1

  ## Set known.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(d.true))    d    = d.true;
  if (!is.null(lambda.true)) lambda = lambda.true;
  if (!is.null(r.true))    r    = r.true;

  ## Preprocess
  ymax = max(y);
  F = cumsum(hist(y, breaks=0:(ymax+1)-0.5, plot=FALSE)$counts)
  G = N - F;
  
  out <- list(beta = matrix(nrow=samp, ncol=P),
                 d = rep(0, samp),
                 lambda = matrix(nrow=samp, ncol=N),
                 r = matrix(nrow=samp, ncol=N)
                 )

  start.time = proc.time()
  
  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if(j==burn+1) start.ess = proc.time()

    ## WARNING: JOINT DRAW.
    ## draw (d, lambda, r | beta)
    ## draw (d | beta)
    phi = drop(X%*%beta)
    mu  = exp(phi)
    d = draw.df(d, mu, G, ymax);
    ## draw (lambda | d, beta)
    psi = phi - log(d);
    p = 1 / (1 + exp(-psi))
    lambda = rgamma(N, y+d, scale=p)
    ## draw (r | d, lambda, beta)
    nmix = compute.mixture(d);
    res  = psi - log(lambda)
    r    = draw.indicators.C(res, nmix);

    ## draw beta
    beta = draw.beta(X, lambda, d, r, nmix, b.0=b.0, P.0=P.0);
    
    # Record if we are past burn-in.
    if (j>burn) {
        out$r[j-burn,]     = r
        out$beta[j-burn,]  = beta
        out$d[j-burn]      = d
        out$lambda[j-burn,]= lambda
    }
    
    if (j %% verbose == 0) { print(paste("NBFS-logmean: Iteration", j)); }
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess

  out
} ## logit.PG.gibbs

################################################################################

if (FALSE) {

  dyn.load("../C/FSF_nmix.so")
  
  N = 500
  X = cbind(1, rnorm(N))

  beta = c(2.0, 1.5);

  ## Poisson-Gamma mixture.
  phi = X %*% beta;
  mu = exp(phi);
  d  = 4
  nmix = compute.mixture(d)
  r = sample.int(nmix$nc, N, replace=TRUE, prob=nmix$p);
  lambda = exp(phi - log(d) - rnorm(N, nmix$m[r], sqrt(nmix$v[r])));
  y = rpois(N, lambda)

  z = log(lambda) + log(d) + nmix$m[r];
  z.hat = z / sqrt(nmix$v[r])
  X.hat = X / sqrt(nmix$v[r])
  lm1 = lm(z.hat ~ X.hat + 0)
  
  samp = 1000
  burn = 500
  verbose = 500
  
  out <- logit.FS.gibbs(y, X,
                        samp=samp, burn=burn, verbose=verbose,
                        beta.true = NULL, d.true=NULL, lambda.true=NULL, r.true=NULL)
  
}


