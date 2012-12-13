
## Independent Metropolis-Hastings and Symmetric Random Walk Metropolis-Hastings.

norm.llh <- function(beta, m, V)
{
  ## beta : P x J-1; coefficients, column assumed to be beta_J = 0.
  ## m : P*(J-1); vectorized mean.
  ## V : P*(J-1) x P*(J-1); Variance matrix for vectorized mean.
  ## P : Precision, in case you do not want to use variance.
  beta = as.numeric(beta)
  ev   = eigen(V);

  P = ev$vectors %*% ( t(ev$vectors) / ev$values );  
  llh = -0.5 * sum( log(abs(ev$values)) ) - 0.5 * t(beta - m) %*% P %*% (beta - m);

  llh
}

norm.llh.2 <- function(beta, m, evalues, P)
{
  ## beta : P x J-1; coefficients, column assumed to be beta_J = 0.
  ## m : P*(J-1); vectorized mean.
  ## V : P*(J-1) x P*(J-1); Variance matrix for vectorized mean.
  ## P : Precision, in case you do not want to use variance.
  ## evalues: evalues of P^{-1}.
  beta = as.numeric(beta)
  
  llh = -0.5 * sum( log(abs(evalues)) ) - 0.5 * t(beta - m) %*% P %*% (beta - m);

  llh
}

std.norm.llh <- function(beta)
{
  beta = as.numeric(beta);
  -0.5 * (beta %*% beta);
}

log.std.norm.dens <- function(ep)
{
  beta = as.numeric(beta);
  d = length(beta);
  -0.5 * d * log(2 * pi) - 0.5 * (ep %*% ep);
}

std.t.llh <- function(beta, df)
{
  ## Assume df is a FIXED parameter.
  beta = as.numeric(beta);
  p = length(beta);
  -0.5 * (df + p) * log(1 + beta %*% beta / df);
}

t.llh <- function(beta, m, V, nu)
{
  beta = as.numeric(beta)
  ev = eigen(V);
  p = length(beta);

  P = ev$vectors %*% ( t(ev$vectors) / ev$values ); 

  llh = lgamma(0.5*(nu+p)) - lgamma(0.5*nu) - 0.5*p*log(nu) - 0.5*p*log(pi) - 
    0.5 * sum( log(abs(ev$values)) ) - 0.5*(nu+p) * log( 1 + t(beta-m) %*% P %*% (beta-m) / nu );

  llh
}

t.llh.2 <- function(beta, m, evalues, P, nu)
{
  beta = as.numeric(beta)
  p = length(beta);

  llh = lgamma(0.5*(nu+p)) - lgamma(0.5*nu) - 0.5*p*log(nu) - 0.5*p*log(pi) - 
    0.5 * sum( log(abs(evalues)) ) - 0.5*(nu+p) * log( 1 + t(beta-m) %*% P %*% (beta-m) / nu);

  llh
}

mlogit.llh <- function(beta, y, X, n,
                       m.0=array(0, dim=c(ncol(X), ncol(y))),
                       P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))))
{
  ## beta : P x J-1; coefficients, column assumed to be beta_J = 0.
  ## y : N x J-1; fraction of responses.
  ## X : N x P; design matrix
  ## P.0 : P x P x J-1; array of matrices for independent prior. 
  beta = matrix(beta, ncol(X), ncol(y));
  XB = X %*% beta;
  ll = sum( n * (rowSums(y * XB) - log( rowSums(exp(XB))+1 )) )
  lr = 0; for (j in 1:ncol(y))
    lr = lr - 0.5 * t(beta[,j] - m.0[,j]) %*% P.0[,,j] %*% (beta[,j] - m.0[,j]);
  ll + lr
}

grad.mlogit.llh <- function(beta, y, X, n,
                       m.0=array(0, dim=c(ncol(X), ncol(y))),
                       P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))))
{
  beta  = matrix(beta, ncol(X), ncol(y));
  expXB = exp(X %*% beta);
  grad  = t(X) %*% ( y - expXB * n / (rowSums(expXB)+1) );
  grad.lr = matrix(0, ncol(X), ncol(y)); for(j in 1:ncol(y))
    grad.lr[,j] = -P.0[,,j] %*% (beta[,j] - m.0[,j]);
  grad = grad + grad.lr
    
  as.numeric(grad)
}

hessian.mlogit.llh <- function(beta, y, X, n, P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))))
{
  beta  = matrix(beta, ncol(X), ncol(y));
  P = ncol(X);
  J = ncol(y) + 1;
  hess = array(0, dim=c(P, J-1, P, J-1));

  XB = X %*% beta
  expXB = exp(XB)
  sum.expXB = rowSums(expXB) + 1;

  for (p in 1:(J-1)) {
    for (j in 1:(J-1)) {

      hess[,j,,p] = t(X) %*% (X * expXB[,p] * expXB[,j] * n / sum.expXB^2)

      if (p == j) {
        hess[,j,,p] = hess[,j,,p] - t(X) %*% (X * expXB[,j] * n / sum.expXB) - P.0[,,j]
      }
            
    }
  }

  hess.vec = array(hess, dim=c(P*(J-1), P*(J-1)));
}

ind.metroplis <- function(y, X, n, m0, P0, m, V, samp=1000, burn=100, tune = 0.25, verbose=1000, df=6)
{
  ## y: response
  ## X: design
  ## n: number of trials per draw.
  ## m: ppsl mean
  ## V: ppsl scale
  
  ## Set Metropolis ##
  evd = eigen(V);
  evalues = evd$values;
  Prec = evd$vectors %*% ( t(evd$vectors) / evd$values );  
  L  = t(chol(V));
  df = abs(df)
  use.t = df!=Inf

  N = nrow(X)
  P = ncol(X)
  J = ncol(y) + 1
  PJ1 = P * (J - 1);

  out <- list(beta = array(0, dim=c(samp, P*(J-1))),
              alpha = rep(0, samp)
              )

  ## start at ppsl mean.
  ep   = rep(0, PJ1)
  beta = m

  ## Set mh.diff
  if (use.t) ppsl.llh = std.t.llh(ep, df) else ppsl.llh = std.norm.llh(ep);
  mh.diff = mlogit.llh(beta, y, X, n, m0, P0) - ppsl.llh

  ## Timing
  start.time = proc.time()
  naccept = 0
  
  ## Do Metropolis ##
  for (i in 1:(samp+burn)) {

    ## beta.0 = beta
    ## optim.out = optim(beta.0, mlogit.llh, gr=grad.mlogit.llh, method="BFGS", hessian=TRUE,
    ##  y=y, X=X, n=n, m.0=m0, P.0=P0, control=list(fnscale=-1));
    
    if (i==burn+1) { start.ess = proc.time(); naccept = 0; }
    
    ## Proposal
    if (use.t) { ep.ppsl = rt(PJ1, df); } else { ep.ppsl = rnorm(PJ1); }
    ppsl = m + tune * (L %*% ep.ppsl);

    ## Ratio -- this works
    ## llh.ratio = mlogit.llh(ppsl, y, X, n, m0, P0) - mlogit.llh(beta, y, X, n, m0, P0);
    ## if (use.t) { prop.ratio = std.t.llh(ep.ppsl, df) - std.t.llh(ep, df); } else
    ## { prop.ratio = std.norm.llh(ep.ppsl) - std.norm.llh(ep); }
    ## log.ratio = llh.ratio - prop.ratio;

    if (use.t) ppsl.llh = std.t.llh(ep.ppsl, df) else ppsl.llh = std.norm.llh(ep.ppsl);
    mh.diff.ppsl = mlogit.llh(ppsl, y, X, n, m0, P0) - ppsl.llh
    log.ratio = mh.diff.ppsl - mh.diff

    ## Accept/Reject
    alpha = min(exp(log.ratio), 1);
    if (runif(1) < alpha) {
      naccept = naccept + 1
      beta = ppsl
      ep   = ep.ppsl
      mh.diff = mh.diff.ppsl
    }
    
    if (i > burn) {
      out$beta[i-burn,] = beta;
      out$alpha[i-burn] = alpha;
    }

    if (i %% verbose == 0) {
      if (i > burn) cat("Ind MH: Ave alpha:", mean(out$alpha[1:(i-burn)]), ", ");
      cat("Iteration:", i, "\n");
    }
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$acceptr    = naccept / samp;

  out
}

sym.rw.metroplis <- function(y, X, n, m0, P0, beta.0, V, samp=1000, burn=100, tune = 0.25, verbose=1000, df=6)
{
  ## y: response
  ## X: design
  ## n: number of trials per draw
  
  N = nrow(X)
  P = ncol(X)
  J = ncol(y) + 1
  PJ1 = P * (J - 1);
  df = abs(df)
  use.t = df!=Inf

  ## Scale of each variable as determined by Hessian.
  sc = sqrt(diag(V));

  out = list(
    beta = array(0, dim=c(samp, P*(J-1))),
    alpha = rep(0, samp)
    )

  beta = beta.0

  ## Timing
  start.time = proc.time()
  
  ## Do Metropolis ##
  for (i in 1:(samp+burn)) {

    if (i==burn+1) start.ess = proc.time();
    
    ## ## Just making stuff up.
    ## H = hessian.mlogit.llh(beta, y, X, n, P0)
    ## evd = eigen(-1*H);
    ## evalues = evd$values;
    ## pos.min = min(evalues[evalues>0]);
    ## evalues[evalues<=0] = pos.min;
    ## RtV = evd$vectors %*% ( t(evd$vectors) / evd$values^0.5 );
    ##  RtV %*%

    ## Proposal -- Symmetric
    if (use.t) { ep = rt(PJ1, df); } else { ep = rnorm(PJ1); }
    ppsl = beta + tune * sc * ep;

    ## Ratio
    log.ratio = mlogit.llh(ppsl, y, X, n, m0, P0) - mlogit.llh(beta, y, X, n, m0, P0)

    ## Accept/Reject?
    alpha = min(exp(log.ratio), 1);
    if (runif(1) < alpha) {
      beta = ppsl
    }
    
    if (i > burn) {
      out$beta[i-burn,] = beta;
      out$alpha[i-burn] = alpha
    }

    if (i %% verbose == 0) {
      if (i > burn) cat("RW MH: Ave alpha:", mean(out$alpha[1:(i-burn)]), ", ");
      cat("Iteration:", i, "\n");
    }
  }
  
  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess

  out
}

mlogit.laplace <- function(y, X, n, m.0, P.0, beta.0=NULL)
{
  y = as.matrix(y)
  X = as.matrix(X)
  
  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;

  if (is.null(beta.0)) beta.0 = matrix(0, P, J-1);
  
  optim.out = optim(beta.0, mlogit.llh, gr=grad.mlogit.llh, method="BFGS", hessian=TRUE,
    y=y, X=X, n=n, m.0=m.0, P.0=P.0, control=list(fnscale=-1));
  
  beta.pm = matrix(optim.out$par, P, J-1);

  ## num - numerically
  ## pen - analytically (pen and paper)
  
  ## hess.num = optim.out$hessian
  hess.pen = hessian.mlogit.llh(beta.pm, y, X, n, P.0)
  
  ## V.num = solve(-1*optim.out$hessian); ## I encountered some numerical instability in the German dataset.
  V.pen = solve(-1*hess.pen);
  ## Chol.ml = chol(V.pen);
  
  m = as.numeric(beta.pm);
  V = V.pen;

  out = list("m"=m, "V"=V)
}

mlogit.MH.R <- function(y, X, n, m.0=NULL, P.0=NULL, beta.0=NULL, samp=1000, burn=1000,
                        method=c("Ind", "RW"), tune=1.0, df=Inf, verbose=1000)
{
  ## psi = X beta where beta = [beta_j], j=1,...,J-1.
  ## J categories, n draws (scalar) for each observation.
  
  ## y: dim N x J-1 response vector of multinomial draw,
  ##    Assume y_i = (y_{i1}, ..., y_{i,J-1}, y_{i,J} = n - ...).
  ## X: dim N x P design matrix -- does not vary with beta_j
  ## m.0: P x J-1 prior mean
  ## P.0: P x P x J-1 prior var.
  ## beta.0: initial value for numerically calculating map.
  
  y = as.matrix(y)
  X = as.matrix(X)
  
  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;
  
  ## n = rep(1, N)
  
  if (is.null(beta.0)) beta.0 = matrix(0, P, J-1);

  if (is.null(m.0)) m.0 = array(0, dim=c(P, J-1));
  if (is.null(P.0)) {
    p.0 = 0e-4
    P.0 = array(diag(p.0, P), dim=c(P, P, J-1));
  }
  
  optim.out = optim(beta.0, mlogit.llh, gr=grad.mlogit.llh, method="BFGS", hessian=TRUE,
    y=y, X=X, n=n, m.0=m.0, P.0=P.0, control=list(fnscale=-1));
  
  beta.pm = matrix(optim.out$par, P, J-1);

  ## num - numerically
  ## pen - analytically (pen and paper)
  
  hess.num = optim.out$hessian
  hess.pen = hessian.mlogit.llh(beta.pm, y, X, n, P.0)
  
  V.num = solve(-1*optim.out$hessian); ## I encountered some numerical instability in the German dataset.
  V.pen = solve(-1*hess.pen);
  ## Chol.ml = chol(V.pen);
  
  m = as.numeric(beta.pm);
  V = V.pen;
  
  if (method[1]=="Ind")
    mh = ind.metroplis(y, X, n, m.0, P.0, m, V, samp=samp, burn=burn, tune=tune, verbose=verbose, df=df)
  if (method[1]=="RW" )
    mh = sym.rw.metroplis(y, X, n, m.0, P.0, beta.pm, V, samp=samp, burn=burn, tune=tune, verbose=verbose, df=df);

  mh$map  = beta.pm;
  mh$var  = V.pen;
  mh$hess = hess.pen;
  
  mh
  
}

################################################################################
                               ## EMPIRICAL KL ##
################################################################################

emp.mlogit.kl <- function(beta, y, X, n, m.0, P.0)
{
  print("THIS DOESN'T WORK.");

  lap.approx = mlogit.laplace(y, X, n, m.0, P.0);
  m = lap.approx$m;
  V = lap.approx$V;
  dim.beta = length(m);
  
  log.mlogit.post <- function(.beta.) { mlogit.llh(.beta., y, X, n, m.0, P.0) }
  mlogit.post <- function(.beta.) { exp(mlogit.llh(.beta., y, X, n, m.0, P.0)) }
  P.log.prop  = mean( apply(beta, 1, log.mlogit.post) )
  log.P.cnst  = log( mean( apply(beta, 1, mlogit.post) ) )

  U = chol(V);
  UInv = backsolve(U, diag(1.0, dim.beta));
  ep = apply(beta, 2, function(.x.) { .x. - m } );
  ep = ep %*% UInv;

  log.norm.post <- function(.ep.) { log.std.norm.dens(.ep.); }
  norm.post <- function(.ep.) { exp(log.std.norm.dens(.ep.)); }
  Q.log.prop = mean( apply(ep, 1, log.norm.post) );
  log.Q.cnst = 0.0; ## log.norm.post is normalized

  KL = P.log.prop - log.P.cnst - Q.log.prop + log.Q.cnst

  out = list("KL"=KL, "P.log.prop"=P.log.prop, "log.P.cnst"=log.P.cnst,
    "Q.log.prop"=Q.log.prop, "log.Q.cnst"=log.Q.cnst);

  out
}

################################################################################
                                   ## MAIN ##
################################################################################

if (FALSE) {

  N = 300;
  P = 2;

  ##------------------------------------------------------------------------------
  ## Correlated predictors
  rho = 0.0
  Sig = matrix(rho, nrow=P, ncol=P); diag(Sig) = 1.0;
  U   = chol(Sig);
  X   = matrix(rnorm(N*P), nrow=N, ncol=P) %*% U;

  beta = rnorm(P, mean=0, sd=2);
  
  psi = X %*% beta;
  p = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);
  n = rep(1, N);
  
}

if (FALSE) {

  y = as.matrix(y)
  X = as.matrix(X)
  
  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;

  n = rep(1, N)
  
  beta.0 = matrix(0, P, J-1);
  
  m.0 = array(0, dim=c(P, J-1));
  p.0 = 0e-4
  P.0 = array(diag(p.0, P), dim=c(P, P, J-1));
  
  optim.out = optim(beta.0, mlogit.llh, gr=grad.mlogit.llh, method="BFGS", hessian=TRUE,
    y=y, X=X, n=n, m.0=m.0, P.0=P.0, control=list(fnscale=-1));
  
  beta.pm = matrix(optim.out$par, P, J-1);
  
  hess.num = optim.out$hessian
  hess.pen = hessian.mlogit.llh(beta.pm, y, X, n, P.0)
  
  V.num = solve(-1*optim.out$hessian);
  V.pen = solve(-1*hess.pen);
  Chol.ml = chol(V.pen);
  
  m = as.numeric(beta.pm);
  V = V.pen;

  samp = 10000
  burn = 1000
  df = Inf
  
  mh.ind = ind.metroplis(y, X, n, m.0, P.0, m, V, samp=samp, burn=burn, tune=1.0, df=df)
  mh.rw  = sym.rw.metroplis(y, X, n, m.0, P.0, beta.pm, V, samp=samp, burn=burn, tune=0.1, df=df);
  ## beta.met.ave = colMeans(mh$beta);
  ## beta.met.ave = array(beta.met.ave, dim=c(P, J-1));
  
  apply(mh.ind$beta[seq(1,samp,2),], 2, sd)
  apply(mh.rw$beta[seq(1,samp,2),] , 2, sd)

  mh = mlogit.MH.R(y, X, n, m.0, P.0, beta.0=beta.pm, method="Ind", tune=1.0, df=Inf)

  emp.mlogit.kl(mh$beta, y, X, n, m.0, P.0);
  
}

################################################################################
                                 ## APPENDIX ##
################################################################################

## For the German Credit Data set it appears that the Hessian is NOT positive
## definite unless I put a very strong prior on it (e.g. p.0 = 1e5).  I'm not
## sure how the algorithm decides to stop when that is the case.

## So why does it work using glm, BayesLogit, or anything else?

## With an improper prior, the Austrailia Credit dataset also not positive
## definite.  It is numerically singular.
