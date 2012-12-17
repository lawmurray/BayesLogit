## Independent Metropolis-Hastings and Symmetric Random Walk Metropolis-Hastings.

################################################################################
      ## These aren't normalized because Ind MH doesn't change the Var ##
################################################################################

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

################################################################################
                                  ## MLOGIT ##
################################################################################

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

mlogit.map <- function(y, X, n, m.0, P.0, beta.0=NULL, calc.V=TRUE)
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
  ## V.pen = solve(-1 * hess.pen);
  ## Chol.ml = chol(V.pen);
  ## V = V.pen;
  
  ## m = as.numeric(beta.pm);

  V = NA;
  if (calc.V) V = solve(-1 * hess.pen);

  out = list("m"=beta.pm, "H"=hess.pen, "V"=V)
}

approx.mlogit <- function(beta, y, X, n, m.0, P.0)
{
  ## Approximate mlogit.post at beta using: - 0.5 h' P h + kapp' h
  grad = grad.mlogit.llh(beta, y, X, n, m.0, P.0)
  Prec = -1 * hessian.mlogit.llh(beta, y, X, n, P.0)
  U = chol(Prec);
  kapp = Prec %*% beta + grad
  m = backsolve(U, kapp, transpose=TRUE)
  m = backsolve(U, m);
  out = list("m"=m, "P"=P, "U"=U, "kapp"=kapp);
  out
}

################################################################################
                                   ## MH 1 ##
################################################################################

ind.metropolis <- function(y, X, n, m0, P0, m, V, samp=1000, burn=100, tune = 0.25, verbose=1000, df=6)
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

  y = as.matrix(y)
  X = as.matrix(X)
  
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

sym.rw.metropolis <- function(y, X, n, m0, P0, beta.0, V, samp=1000, burn=100, tune = 0.25, verbose=1000, df=6)
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
    mh = ind.metropolis(y, X, n, m.0, P.0, m, V, samp=samp, burn=burn, tune=tune, verbose=verbose, df=df)
  if (method[1]=="RW" )
    mh = sym.rw.metropolis(y, X, n, m.0, P.0, beta.pm, V, samp=samp, burn=burn, tune=tune, verbose=verbose, df=df);

  mh$map  = beta.pm;
  mh$var  = V.pen;
  mh$hess = hess.pen;
  
  mh
  
}

################################################################################
                                 ## MH 2 ##
################################################################################

llh.norm <- function(x, kappa, UorP, is.prec=TRUE, is.mean=TRUE)
{
  ## x: point at which llh is evaluated.
  ## kappa: Prec * mean or kappa = mean.
  ## UorP: chol(Prec) or Prec
  ## is.prec: specifies that UorP is Prec.
  ## is.mean: kappa is mean
  
  if (is.prec) U = chol(UorP)
  else U = UorP;

  if (is.mean) m = kappa
  else {
    m = backsolve(U, kappa, transpose=TRUE);
    m = backsolve(U, m);
  }
  e = x - m;
  e = U %*% e;

  -1 * sum( log(diag(U)) ) - 0.5 * t(e) %*% e;
}

llh.norm.utest <- function(x1=c(1,2), x2=c(-1.32,1), m=c(0.5,0.5), V=diag(1, 2))
{  
  m = as.matrix(m)
  V = as.matrix(V)
  
  require("mvtnorm")
  k = V %*% m;
  P = solve(V);

  d1 = llh.norm(x1, k, P) - llh.norm(x2, k, P);
  d2 = dmvnorm(x1, m, V, log=TRUE) - dmvnorm(x2, m, V, log=TRUE)

  c(d1,d2)
}

llh.t <- function(x, kappa, UorP, df, is.prec=TRUE, is.mean=TRUE)
{
  ## Assume the df are FIXED between comparisons.
  ## x: point at which llh is evaluated.
  ## kappa: Prec * mean or kappa = mean.
  ## UorP: chol(Prec) or Prec
  ## df: degrees of freedom
  ## is.prec: specifies that UorP is Prec.
  
  if (is.prec) U = chol(UorP)
  else U = UorP;

  p = length(x);

  if (is.mean) m = kappa
  else {
    m = backsolve(U, kappa, transpose=TRUE);
    m = backsolve(U, m);
  }
  e = x - m;
  e = U %*% e;

  -1 * sum( log(diag(U)) ) - 0.5 * (df + p) * log(1 + t(e) %*% e / df);
}

llh.t.utest <- function(x1=c(1,2), x2=c(-1.32,1), m=c(0.5,0.5), V=diag(1, 2), df=6)
{  
  m = as.matrix(m)
  V = as.matrix(V)
  
  require("mvtnorm")
  k = V %*% m;
  P = solve(V);

  d1 = llh.t(x1, k, P, df) - llh.t(x2, k, P, df);
  d2 = dmvt(x1, m, V, df=df, log=TRUE) - dmvt(x2, m, V, df=df, log=TRUE)

  c(d1,d2)
}

log.tget.to.ppsl <- function(beta, y, X, n, m.0, P.0, m, UorP, df=Inf, is.prec=TRUE)
{
  ## Calculate log f - log q
  if (df==Inf)
    log.fdivq = llh.norm(beta, m, UorP, is.prec)
  else
    log.fdivq = llh.t(beta, m, UorP, df, is.prec)
  
  log.fdivq = mlogit.llh(beta, y, X, n, m.0, P.0) - log.fdivq

  log.fdivq
}

draw.beta.ind.MH <- function(beta, y, X, n, m0, P0, map=NULL, U.map=NULL, df=Inf, log.fdivq=NULL)
{
  ## Ind MH sample when m0, P0 might be changing.
  ## beta: previous beta
  ## y, X, n: data
  ## m0, P0: priors for beta
  ## map: posterior mode
  ## U.map: chol(-H) where H is Hessian at mode
  ## df: degrees of freedom for proposal
  ## log.fdiq.q : previous
  
  P = length(beta)
  zero = rep(0, P);

  ## Calculate Laplace Approx if necessary.
  if (recalc <- (is.null(map) || is.null(U.map))) {
    lap = mlogit.map(y, X, n, m.0, P.0, beta.0=NULL, calc.V=FALSE)
    map   = lap$m
    P.map = -1 * lap$H
    U.map = chol(P.map)
  }
  
  ## Calculate old log f - log q
  if (is.null(log.fdivq) || recalc) {
    if (df==Inf)
      log.fdivq = llh.norm(beta, map, U.map, is.prec=FALSE)
    else
      log.fdivq = llh.t(beta, map, U.map, df, is.prec=FALSE)

    log.fdivq = mlogit.llh(beta, y, X, n, m.0, P.0) - log.fdivq
  }

  ## log.fdivq = log.tget.to.ppsl(beta, y, X, n, m.0, P.0, map, U.map, df=df, is.prec=FALSE)

  ## Propose
  if (df==Inf) ep.ppsl = rnorm(P) else ep.ppsl = rt(P, df=df)
  ppsl = map + backsolve(U.map, ep.ppsl);

  ## Claculate new log f - log q
  if (df==Inf)
    log.fdivq.ppsl = llh.norm(ppsl, map, U.map, is.prec=FALSE)
  else
    log.fdivq.ppsl = llh.t(ppsl, map, U.map, df, is.prec=FALSE)
  
  log.fdivq.ppsl = mlogit.llh(ppsl, y, X, n, m.0, P.0) - log.fdivq.ppsl

  ## acceptance prob
  alpha = min( exp(log.fdivq.ppsl - log.fdivq), 1);
  accept = runif(1) < alpha
  
  if (accept) {
    beta = ppsl;
    log.fdivq = log.fdivq.ppsl
  }

  out = list("beta"=beta, "alpha"=alpha, "accept"=accept, "log.fdivq"=log.fdivq)

  out
}

draw.beta.rw.MH <- function(beta, y, X, n, m.0, P.0, df=Inf)
{
  ## Make draw based upon Taylor Approximation at beta
  
  y = as.matrix(y)
  X = as.matrix(X)
  
  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;

  ## Calculate Taylor Approx at beta.
  approx = approx.mlogit(beta, y, X, n, m.0, P.0)
  m  = approx$m
  U  = approx$U

  ## Propose
  if (df==Inf) ep.ppsl = rnorm(P) else ep.ppsl = rt(P, df=df)
  ppsl = m + backsolve(U, ep.ppsl);

  ## log.fdivq.ppsl
  log.fdivq.ppsl = log.tget.to.ppsl(ppsl, y, X, n, m.0, P.0, m, U, df=df, is.prec=FALSE)

  ## Calculate Taylor Approx at ppsl.
  approx = approx.mlogit(ppsl, y, X, n, m.0, P.0)
  m  = approx$m
  U  = approx$U
  
  ## log.fdivq.beta
  log.fdivq.beta = log.tget.to.ppsl(beta, y, X, n, m.0, P.0, m, U, df=df, is.prec=FALSE)

  ## acceptance prob
  alpha = min( exp(log.fdivq.ppsl - log.fdivq.beta), 1);
  accept = runif(1) < alpha;
  
  if (accept) {
    beta = ppsl;
  }

  out = list("beta"=beta, "alpha"=alpha, "accept"=accept)

  out
}

ind.metropolis.2 <- function(y, X, n, m.0, P.0, samp=1000, burn=100, verbose=1000, df=6)
{
  ## An independent MH routine that can be adapted to varying P0.
  ## Assume proper posterior
  
  ## y: response, number of responses in each category
  ## X: design
  ## n: number of trials per draw.
  
  y = as.matrix(y)
  X = as.matrix(X)
  
  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;
  PJ1 = P * (J - 1);
  
  ## n = rep(1, N)
  
  if (is.null(m.0)) m.0 = array(0, dim=c(P, J-1));
  if (is.null(P.0)) {
    p.0 = 0e-4
    P.0 = array(diag(p.0, P), dim=c(P, P, J-1));
  }

  out <- list(beta = array(0, dim=c(samp, P*(J-1))),
              alpha = rep(0, samp)
              )

  ## Find mode
  out.map = mlogit.map(y, X, n, m.0, P.0, beta.0=NULL)
  m = out.map$m
  U = chol(-1 * out.map$H)
  
  ## start at ppsl mean.
  ep   = rep(0, PJ1)
  beta = out.map$m

  ## Generate proposal
  out.beta = draw.beta.ind.MH(beta, y, X, n, m.0, P.0, map=m, U.map=U, df=df, log.fdivq=NULL)
  beta = out.beta$beta
  log.fdivq = out.beta$log.fdivq

  ## Timing
  start.time = proc.time()
  naccept = 0
  
  ## Do Metropolis ##
  for (i in 1:(samp+burn)) {

    if (i==burn+1) { start.ess = proc.time(); }

    out.beta = draw.beta.ind.MH(beta, y, X, n, m.0, P.0, map=m, U.map=U, df=df, log.fdivq=log.fdivq)
    ## out.beta = draw.beta.ind.MH(beta, y, X, n, m.0, P.0, map=NULL, U.map=NULL, df=df, log.fdivq=NULL)
    beta = out.beta$beta
    log.fdivq = out.beta$log.fdivq
    
    if (i > burn) {
      out$beta[i-burn,] = beta
      out$alpha[i-burn] = out.beta$alpha
      naccept = naccept + out.beta$accept
    }

    if (i %% verbose == 0) {
      if (i > burn) cat("Ind MH2: Ave alpha:", mean(out$alpha[1:(i-burn)]), ", ");
      cat("Iteration:", i, "\n");
    }
    
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$acceptr    = naccept / samp;

  out
}

rw.metropolis <- function(y, X, n, m.0=NULL, P.0=NULL, samp=1000, burn=100, verbose=1000, df=6)
{
  ## RW MH based upon Taylor approx. to post. at current beta.
  ## Assume proper posterior
  
  ## y: response, number of responses in each category
  ## X: design
  ## n: number of trials per draw.
  
  y = as.matrix(y)
  X = as.matrix(X)
  
  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;
  PJ1 = P * (J - 1);
  
  ## n = rep(1, N)
  
  if (is.null(m.0)) m.0 = array(0, dim=c(P, J-1));
  if (is.null(P.0)) {
    p.0 = 0e-4
    P.0 = array(diag(p.0, P), dim=c(P, P, J-1));
  }

  out <- list(beta = array(0, dim=c(samp, P*(J-1))),
              alpha = rep(0, samp)
              )

  ## Find mode
  out.map = mlogit.map(y, X, n, m.0, P.0, beta.0=NULL)
  m = out.map$m
  U = chol(-1 * out.map$H)
  
  ## start at ppsl mean.
  ep   = rep(0, PJ1)
  beta = out.map$m

  ## Generate proposal
  out.beta = draw.beta.rw.MH(beta, y, X, n, m.0, P.0, df=df)
  beta = out.beta$beta

  ## Timing
  start.time = proc.time()
  naccept = 0
  
  ## Do Metropolis ##
  for (i in 1:(samp+burn)) {

    if (i==burn+1) { start.ess = proc.time(); }

    out.beta = draw.beta.rw.MH(beta, y, X, n, m.0, P.0, df=df)
    beta = out.beta$beta
    
    if (i > burn) {
      out$beta[i-burn,] = beta
      out$alpha[i-burn] = out.beta$alpha
      naccept = naccept + out.beta$accept
    }

    if (i %% verbose == 0) {
      if (i > burn) cat("RW MH2: Ave alpha:", mean(out$alpha[1:(i-burn)]), ", ");
      cat("Iteration:", i, "\n");
    }
    
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$acceptr    = naccept / samp;

  out
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

  N = 100;
  P = 2;

  ##------------------------------------------------------------------------------
  ## Correlated predictors
  rho = 0.99
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
  
  mh.ind = ind.metropolis(y, X, n, m.0, P.0, m, V, samp=samp, burn=burn, tune=1.0, df=df)
  mh.rw  = sym.rw.metropolis(y, X, n, m.0, P.0, beta.pm, V, samp=samp, burn=burn, tune=0.1, df=df);
  mh.ind2 = ind.metropolis.2(y, X, n, m.0, P.0, samp=samp, burn=burn, df=df);
  mh.rw2  = rw.metropolis(y, X, n, m.0, P.0, samp=samp, burn=burn, df=df);
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

## It appears for the second MH sampler that things are much faster without
## function calls.  I think this makes sense because R is always copying objects
## when passing to a function.
