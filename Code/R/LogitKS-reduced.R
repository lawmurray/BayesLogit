## This follows Holmes and Held (2006).

################################################################################

if (!is.loaded("hh_lambda")) dyn.load("~/RPackage/BayesLogit/Code/R/hh_lambda.so")

draw.lambda.C <- function(N, r)
{
  ok = all(r>0);
  if (!ok) {
    print("r must be > 0.")
    return(NULL);
  }
  
  r      = array(r, N);
  lambda = rep(0.0, N)
  iter   = rep(0.0, N);
  
  ## for (i in 1:N) {
  ##   OUT = .C("hh_lambda", lambda[i], r[i], iter[i]);
  ##   lambda[i] = OUT[[1]];
  ##   iter[i]   = OUT[[3]];
  ## }
  ## lambda
  
  OUT = .C("hh_lambda_vec", lambda, r, N); 
  OUT[[1]]
}

################################################################################

draw.beta <- function(X, Z, lambda, P.0=0.0)
{
  N = nrow(X);
  P = ncol(X);

  ## if (is.null(P.0)) P.0 = matrix(0, P, P);

  ## Posterior precision, variance, and mean.
  P.1 = t(X) %*% (X / lambda) + P.0;
  
  V.1 = try(solve(P.1))
  if (!is.matrix(V.1)) {
    print(eigen(P.1)$values);
    print(lambda);
    print(P.1[1:10,1:10]);
  }
  m.1 = V.1 %*% (t(X) %*% (Z / lambda));

  beta = m.1 + t(chol(V.1)) %*% rnorm(P);
} ## draw.beta

## Easy to draw logistic by inversion.
rtlogis.upper <- function(N, mu, s, upper)
{
  V = runif(N) * plogis(upper, mu, s);
  X = mu - s * log(1/V - 1)
}

rtlogis.lower <- function(N, mu, s, lower)
{
  p = plogis(lower, mu, s, lower.tail=FALSE);
  V = (1-p) + runif(N) * p
  X = mu - s * log(1/V - 1)
}

draw.z.lambda <- function(y, Xbeta)
{
  N = length(y);

  ## Draw Z.
  Z = rep(0.0, N);

  Z1 = rtlogis.lower(N, Xbeta, 1.0, lower=0.0);
  Z2 = rtlogis.upper(N, Xbeta, 1.0, upper=0.0);
  Z = y * Z1 + (1-y) * Z2;

  ## Draw lambda.
  r = abs(Z-Xbeta);
  lambda = draw.lambda.C(N, r);

  data.frame("Z"=Z, "lambda"=lambda)
} ## draw.z.lambda

draw.Z <- function(y, Xbeta, lambda)
{
  N = length(Xbeta);
  
  Z = rep(0.0, N);

  Z1 = rtnorm(N, Xbeta, sqrt(lambda), lower=0.0)
  Z2 = rtnorm(N, Xbeta, sqrt(lambda), upper=0.0)
  Z  = y * Z1 + (1-y) * Z2;
} ## draw.Z

logit.KS.gibbs <- function(y, X, samp=1000, burn=100, P.0=NULL, verbose=1000,
                           beta.true=NULL, Z.true=NULL, lambda.true=NULL)
{
  N = nrow(X);
  P = ncol(X);
  M = samp;
  
  out = list(
    beta   = array(0, dim=c(M, P)),
    Z      = array(0, dim=c(M, N)),
    lambda = array(0, dim=c(M, N))
    )

  if(is.null(P.0)) P.0 = matrix(0, P, P);
  
  ## Initilaize Sampler.
  lambda = rep(1.0, N);
  W = rlogis(N, 0.0, 1.0);
  Z = W * (W >= 0 & y==1) + W * (W < 0 & y==0);
  beta = draw.beta(X, Z, lambda, P.0)

  ## Set by known
  if (!is.null(beta.true))   beta   = beta.true
  if (!is.null(Z.true))      Z      = Z.true
  if (!is.null(lambda.true)) lambda = lambda.true

  for (i in 1:(samp+burn)) {

    ## BEGIN DRAW BETA
    ## beta = draw.beta(X, Z, lambda, P.0)
    P.1 = t(X) %*% (X / lambda) + P.0;
    V.1 = solve(P.1)
    ## V.1 = try(solve(P.1))
    ## if (!is.matrix(V.1)) {
    ##   print(eigen(P.1)$values);
    ##   print(lambda);
    ##   print(P.1[1:10,1:10]);
    ## }
    m.1 = V.1 %*% (t(X) %*% (Z / lambda));

    beta = m.1 + t(chol(V.1)) %*% rnorm(P);
    ## END DRAW BETA
    
    Xbeta = X %*% beta;


    ## BEGIN DRAW Z, LAMBDA
    ## zlambda = draw.z.lambda(y, Xbeta);
    ## Z = zlambda$Z;
    ## lambda = zlambda$lambda;

    mu = Xbeta;
    s  = 1.0
    p = plogis(0.0, mu, s, lower.tail=FALSE);
    V = (1-p) + runif(N) * p
    Z1 = mu - s * log(1/V - 1)

    mu = Xbeta;
    s  = 1.0;
    V = runif(N) * plogis(0.0, mu, s);
    Z2 = mu - s * log(1/V - 1)

    Z = y * Z1 + (1-y) * Z2;

    ## Draw lambda.
    r = abs(Z-Xbeta);
    lambda = draw.lambda.C(N, r);
    ## END DRAW Z, LAMBDA
    
    if (i > burn) {
      out$beta  [i-burn,] = beta;
      out$Z     [i-burn,] = Z;
      out$lambda[i-burn,] = lambda;
    }

    if (i %% verbose == 0) cat("Iteration", i, "\n");
    
  }

  out
}

################################################################################
                            ## Utility Functions ##
################################################################################

plot.compare <- function(out.hh, out.bl)
{
  P = ncol(out.hh$beta)
  old.par = par(mfrow=c(2,2));
  
  for (i in 1:P) {
    hist(out.hh$beta[,i], main="H&H", xlab=paste("beta", i));
    acf(out.hh$beta[,i], main="H&H");
    cat("mean:", mean(out.hh$beta[,i]), "sd:", sd(out.hh$beta[,i]), "\n");
    hist(out.bl$beta[,i], main="BL", xlab=paste("beta", i));
    acf(out.bl$beta[,i], main="BL");
    cat("mean:", mean(out.bl$beta[,i]), "sd:", sd(out.bl$beta[,i]), "\n");
    readline("Press <ENTER> to continue.");
  }

  par(mfrow=old.par);
}

################################################################################
                                   ## MAIN ##
################################################################################

## Synthetic test 1 ##

if (FALSE) {
  source("KS.R")
  
  N = 1000;
  samp = 1000;
  burn = 500;
  
  beta = c(1.0, 0.4);
  X = cbind(1, rnorm(N));
  psi = rks(N);
  lambda = (2*psi)^2;
  ep = rnorm(N, 0.0, sqrt(lambda));
  z  = X %*% beta + ep
  y  = as.numeric(z > 0);

  start = proc.time()
  out = logit.KS.gibbs(y, X, samp=samp, burn=0, verbose=100,
                       beta.true=beta, Z.true=z, lambda.true=lambda);
  total.time = proc.time() - start;

  glm1   = glm(y ~ X+0, family=binomial(link="logit"));
  summary(glm1)
  ## out.bl = logit.gibbs(y, X, samp=300, burn=0, n.prior=0);

}


################################################################################
                                 ## APPENDIX ##
################################################################################

## Upper Truncation
  ## kern: \1_{(-\infty, u)} f(x) dx.
  ## \int^a kern: F(\min{a,u})
  ## C = F(u)
  ## CDF: F(min(a,u))/F(u)

## Lower Truncation
  ## kern: \1_{(u, \infty)} f(x) dx
  ## \int^a kern: 0, if a < u; F(a) - F(u), else
  ## C = 1 - F(u)
  ## CDF: \1{a > u} [ F(a) - F(u) ] / [ 1 - F(u) ]
