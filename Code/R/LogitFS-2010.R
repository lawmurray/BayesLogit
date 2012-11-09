## This R script implements Fruhwirth-Schnatter and Fruhwirth's normal-mixture
## approximation to logistic regression (2010).

################################################################################

## Based on Monahan & Stefanski means based on their method.  They do not use
## the numbers in their article.  I'm not sure what it is.  They use the same
## probabilities but different variances?

## ## Define normal mixture -- 6 comp.  FS&F p. 119, based on Monahan & Stefanski (1992).
## normal.mixture = list(
##   w = c(1.8446, 17.268, 37.393, 31.697, 10.89, 0.90745) / 100,
##   m = rep(0, 6),
##   v = c(0.68159, 1.2419, 2.2388, 4.0724, 7.4371, 13.772)
##   )
##   normal.mixture$s = sqrt(normal.mixture$v)
##   normal.mixture$N = length(normal.mixture$w)
## c(1.21126, 0.89734, 0.66833, 0.49553, 0.36668, 0.26946)

## ## Define normal mixture -- 3 comp.  FS&F p. 119, based on Monahan & Stefanski (1992).
## normal.mixture = list(
##   w = c(25.22, 58.523, 16.257) / 100,
##   m = rep(0, 3),
##   v = c(1.2131, 2.9955, 7.5458)
##   )
##   normal.mixture$s = sqrt(normal.mixture$v)
##   normal.mixture$N = length(normal.mixture$w)

## Define normal mixture -- 6 comp.  FS&F p. 119, based on K-L distance.
normal.mixture = list(
  w = c(5.8726, 28.74, 36.756, 22.427, 5.8701, 0.33466) / 100,
  m = rep(0, 6),
  v = c(0.84678, 1.61, 2.8904, 5.0772, 8.9109, 15.923)
  )
  normal.mixture$s = sqrt(normal.mixture$v)
  normal.mixture$N = length(normal.mixture$w)

## lgs.samp = -log(rexp(10000)) + log(rexp(10000));

## ## Define normal mixture -- 8 comp.  Monahan & Stefanski (1992).  They use
## ## inverse-scale.  NOTE! Initially, I only used the first 5-digits.  That was
## ## WAY off (in terms of the estimated variance).  You need a lot of precision.
## normal.mixture = list(
##   w = c(0.0032463432, 0.0515174770, 0.1950779126, 0.3155698236,
##         0.2741495761, 0.1310768806, 0.0279241871, 0.0014495678),
##   m  = rep(0, 8),
##   is = c(1.3653408062, 1.0595239710, 0.8307913137, 0.6507321666,
##          0.5081354253, 0.3963122451, 0.3089042522, 0.2382126164)
##   )
##   normal.mixture$s = 1 / normal.mixture$is
##   normal.mixture$v = normal.mixture$s^2
##   normal.mixture$N = length(normal.mixture$w)

lgs.samp = rlogis(10000)
normal.mixture$mar.mean = 0
normal.mixture$mar.var  = pi^2 / 3;

## Make a copy.
NM = normal.mixture

################################################################################

################################################################################

draw.beta <- function(z, X, r, b.0=NULL, B.0=NULL, P.0=NULL)
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
  a.L = t(Xdv) %*% z;

  P.N = P.0 + P.L;
  ## B.N = solve(P.N);
  B.N = chol2inv(chol(P.N));
  b.N = B.N %*% (a.L + P.0 %*% b.0);
  
  beta = b.N + t(chol(B.N)) %*% rnorm(P)
} ## draw.beta

draw.z <- function(lambda, y){
  n = length(lambda)
  u = runif(n)
  z = log(lambda * u + y) - log(1 - u + lambda * (1-y));
  z
} ## draw.z

logit.mix.gibbs <- function(y, X, samp=1000, burn=100, b.0=NULL, B.0=NULL, P.0=NULL, verbose=10000,
                            beta.true=NULL, z.true=NULL, r.true=NULL)
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
    r    = array(0, dim=c(M, N))
    )

  beta = matrix(0.0, P);  # even odds.
  z    = matrix(0.0, N);
  r    = matrix(1, N);

  ## Always be careful with indicators.  You don't want to start with a very
  ## unlikely value (given the data).  It doesn't matter here because we can
  ## sample y.u without respect to r to start.  Nonetheless, seed randomly out
  ## of principle.
  r = sample.int(NM$N, N, prob=NM$w, replace=TRUE);

  ## In case we are doing testing.  May remove later.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(z.true))    z = z.true;
  if (!is.null(r.true))    r = r.true;

  start.time = proc.time()
  
  for (i in 1:(samp+burn)) {
    if (i==burn+1) start.ess = proc.time()

    lambda = drop(exp(X %*% beta));

    ## WARNING: (z | r, beta, y) != (z | beta, y).
    ## JOINT DRAW: (z, r | beta, y) = (z | beta, y) (r | z, beta, y)
    z    = draw.z(lambda, y);
    r    = draw.indicators.logis.C(z, lambda, nmix)
    ## I tried inserting the code directly.  It did not speed things up.
    
    ## (beta | r, z, y)
    beta = draw.beta(z, X, r, b.0, "P.0"=P.0);
    
    if (i > burn) {
      out$beta[i-burn,] = beta;
      out$z [i-burn,] = z;
      out$r   [i-burn,] = r;
    }

    if (i %% verbose == 0) cat("LogitFS-2010: Iteration", i, "\n");  
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess

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
  ep   = rnorm(N, rep(0, N), NM$s[r]);
  z    = X %*% beta + ep
  lambda = exp(X %*% beta);
  y    = as.numeric(z > 0);
  
  B.0 = diag(1, 2)
  out.fs = logit.mix.gibbs(y, X, samp=samp, burn=burn, beta.true=beta, z.true=z, r.true=r);

  glm1   = glm(y ~ X+0, family=binomial(link="logit"));
  print(coef(glm1));
  out.bl = logit(y, X, samp=samp, burn=burn, n.prior=0);

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

################################################################################
                                 ## APPENDIX ##
################################################################################

## I tried drawing z and r within the Gibbs loop, but that didn't save any time.
