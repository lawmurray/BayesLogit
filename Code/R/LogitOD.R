library("BayesBridge")
## library("msm")
library("coda")

logit.OD.R <- function(y, X, m0=NULL, P0=NULL, samp=1000, burn=100, verbose=100,
                       beta.true=NULL, phi.true=NULL, z.true=NULL)
{
  ## y: binary response
  ## X: design matrix
  
  X = as.matrix(X)

  N = nrow(X)
  P = ncol(X)

  ## Storage
  out <- list(z = matrix(0, nrow=samp, ncol=N),
              phi = matrix(0, nrow=samp, ncol=N),
              beta = matrix(0, nrow=samp, ncol=P))

  ## Prior
  if (is.null(m0)) m0 = rep(0, P)
  if (is.null(P0)) P0 = diag(0, P)
  
  ## Constants
  nu = 7.3
  tilsig2 = pi^2 * (nu - 2) / (3 * nu)
  left =rep(0, N); left [y==0]=-Inf;
  right=rep(0, N); right[y==1]=Inf;
  k0 = P0 %*% m0;

  ## Track
  z    = rep(0, N)
  beta = rep(0, P)
  phi  = 1;

  ## Known
  if (!is.null(z.true)) z = z.true
  if (!is.null(phi.true)) phi = phi.true
  if (!is.null(beta.true)) beta = beta.true

  start.time = proc.time()
  start.ess  = proc.time()
  
  for (i in 1:(samp+burn)) {

    if (i==burn) start.ess = proc.time()
    
    psi = X %*% beta;

    ## z
    z = BayesBridge::rtnorm(N, psi, sqrt(tilsig2 / phi), left, right)

    ## phi
    res   = z - psi;
    shape = 0.5 * (nu + 1);
    rate  = 0.5 * (nu + res^2 / tilsig2);
    phi = rgamma(N, shape, rate=rate);

    ## beta
    PP = t(X) %*% (X * phi) / tilsig2 + P0;
    U  = chol(PP)
    PV = chol2inv(U)
    mm = PV %*% (k0 + t(X) %*% (z * phi) / tilsig2);
    beta = mm + t(chol(PV)) %*% rnorm(P)

    if (i > burn) {
      out$beta[i-burn,] = beta
      out$z[i-burn,]    = z
      out$phi[i-burn,]  = phi      
    }

    end.time = proc.time();
    out$total.time = end.time - start.time
    out$ess.time   = end.time - start.ess
    
    if (i %% verbose == 0) cat("LogitOD: Iteration", i, "\n");
    
  }

  out
}

################################################################################

llh.od <- function(beta, y, X, nu, sig2=1)
{
  sig = sqrt(sig2)
  psi = X %*% beta
  p = pt(psi/sig, nu)
  llh = sum(y * log(p) + (1-y) * log(1-p))
  llh
}

od.map <- function(y, X, nu=7, sig2=1, beta.0=NULL)
{
  X = as.matrix(X)
  N = nrow(X)
  P = ncol(X)
  if (is.null(beta.0)) beta.0 = rep(0, P)
  
  optim.out = optim(beta.0, llh.od, method="BFGS", hessian=TRUE,
    y=y, X=X, nu=nu, sig2=sig2, control=list(fnscale=-1));

  m = optim.out$par
  H = optim.out$hessian
  V = solve(-H)

  out = list(H=H, V=V, m=m);

  out
}

################################################################################

if (FALSE) {
  
  N = 5000
  nu = 7.3
  tilsig2 = pi^2 * (nu - 2) / (3 * nu)

  X = cbind(1, rnorm(N));

  beta = c(1.0, 1.5);
  psi = X %*% beta;
  phi = rgamma(N, nu/2, rate=nu/2)
  z   = rnorm(N, psi, sqrt(tilsig2 / phi))
  y   = as.numeric(z > 0);
  
  ## p = exp(psi) / (1 + exp(psi));
  ## y = rbinom(N, 1, prob=p);

  map.out = od.map(y, X, nu, tilsig2)
  
  samp=1000
  burn=200
  verbose=100

  ## source("LogitOD.R")
  out.od <- logit.OD.R(y, X, m0=NULL, P0=NULL, samp=samp, burn=burn, verbose=verbose,
                       beta.true=NULL, phi.true=NULL, z.true=NULL)
  glm1 = glm(formula = y ~ X + 0, family = binomial(link = logit))
  
  colMeans(out.od$beta)
  coef(glm1)
  apply(out.od$beta, 2, effectiveSize)
  
}
