
## binary logit with groups with different intercepts.

################################################################################
                             ## For mixed models ##
################################################################################

blogit.llh.mm <- function(abphi, y, X.re, X.fe, n, shape=1, rate=1,
                          m.0=array(0, dim=c(ncol(X.fe))),
                          P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))))                         
{
  ## abphi: (P.a + P.b + 1) x J-1; coefficients, column assumed to be beta_J = 0.
  ## y : N; number of responses
  ## X : N x P.ab: design matrix
  ## P.0 : P.b x P.b; array of matrices for independent prior. 
  
  N = length(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  X   = cbind(X.re, X.fe)
  P.a = ncol(X.re)
  P.b = ncol(X.fe)
  P.ab = P.a + P.b
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  phi.idc = P.ab + 1;

  alpha = abphi[a.idc]
  beta = abphi[b.idc]
  phi  = abphi[phi.idc]

  if (phi < 0) return(-Inf)
  
  ab = abphi[1:P.ab]
  Psi = X %*% ab;
  ## p = exp(Psi) / (1 + exp(Psi));

  llh = t(y) %*% Psi - sum(n * log(1 + exp(Psi)));
  llh = llh - 0.5 * t(beta - m.0) %*% P.0 %*% (beta - m.0);
  llh = llh + (0.5 * (P.a + shape) - 1) * log(phi) - 0.5 * phi * ( sum(alpha*alpha) + rate );

  llh
}

blogit.llh.mm.alt <- function(abphim, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                              m.0=array(0, dim=c(ncol(X.fe))),
                              P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))))                         
{
  ## abphi: (P.a + P.b + 1) x J-1; coefficients, column assumed to be beta_J = 0.
  ## y : N; number of responses
  ## X : N x P.ab: design matrix
  ## P.0 : P.b x P.b; array of matrices for independent prior. 
  
  N = length(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  X   = cbind(X.re, X.fe)
  P.a = ncol(X.re)
  P.b = ncol(X.fe)
  P.ab = P.a + P.b
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  phi.idc = P.ab + 1;

  alpha = abphim[a.idc]
  beta = abphim[b.idc]
  phi  = abphim[phi.idc]
  m    = abphim[phi.idc+1]

  if (phi < 0) return(-Inf)
  
  ab = abphim[1:P.ab]
  Psi = X %*% ab;
  ## p = exp(Psi) / (1 + exp(Psi));

  llh = t(y) %*% Psi - sum(n * log(1 + exp(Psi)));
  llh = llh - 0.5 * t(beta - m.0) %*% P.0 %*% (beta - m.0);
  llh = llh + (0.5 * (P.a + shape + 1) - 1) * log(phi) - 0.5 * phi * ( sum((alpha-m)^2) + m^2/kappa^2 + rate );

  llh
}

## DO NOT USE!!!
## I think something may be wrong here.
grad.blogit.llh.mm <- function(abphi, y, X.re, X.fe, n, shape=1, rate=1,
                               m.0=array(0, dim=c(ncol(X.fe))),
                               P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))) )
                               
{
  N = length(y)
  X = cbind(X.re, X.fe);
  P.a = ncol(X.re)
  P.b = ncol(X.fe)
  P.ab = P.a + P.b
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  phi.idc = P.ab + 1;

  alpha = abphi[a.idc]
  beta = abphi[b.idc]
  phi  = abphi[phi.idc]
  
  ab = abphi[1:P.ab]
  Psi = X %*% ab;
  p = exp(Psi) / (1 + exp(Psi));

  grad = rep(0, P.ab+1)
  grad[1:P.ab] = t(y - p * n) %*% X;
  grad[a.idc] = grad[a.idc] - phi * alpha;
  grad[phi.idc] = ( 0.5 * (P.a + shape) - 1 ) / phi - 0.5 * (sum(alpha*alpha) + rate);
  
  grad
}

hessian.blogit.llh.mm <- function(abphi, y, X.re, X.fe, n, shape=1, rate=1,
                                  ## m.0=array(0, dim=c(ncol(X.fe))),
                                  P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))) )

{
  N = length(y)
  X = cbind(X.re, X.fe);
  P.a = ncol(X.re)
  P.b = ncol(X.fe)
  P.ab = P.a + P.b
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  phi.idc = P.ab + 1;

  alpha = abphi[a.idc]
  beta = abphi[b.idc]
  phi  = abphi[phi.idc]
  
  ab = abphi[1:P.ab]
  Psi = X %*% ab;
  p = exp(Psi) / (1 + exp(Psi));

  d = (p^2 - p) * n
  H = t(X) %*% (X * rep(d, P.ab));
  
  hess = matrix(0, P.ab+1, P.ab+1);
  hess[1:P.ab, 1:P.ab] = H
  hess[a.idc, a.idc] = hess[a.idc, a.idc] - diag(phi, P.a);
  hess[a.idc, phi.idc] = -1 * alpha;
  hess[phi.idc, a.idc] = -1 * alpha;
  hess[phi.idc, phi.idc] = -1 * (0.5 * (P.a + shape) - 1) / phi^2;

  hess
}

blogit.mm.map <- function(y, X.re, X.fe, n, shape=1, rate=1,
                       m.0=array(0, dim=c(ncol(X.fe))),
                       P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                       abphi.0=NULL, calc.V=TRUE, trace=FALSE)
{ 
  X = cbind(X.re, X.fe);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(abphi.0)) { abphi.0 = matrix(0, P+1); abphi.0[P+1] = 1; }

  optim.out <- optim(abphi.0, blogit.llh.mm, hessian=FALSE,
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate,
                     m.0=m.0, P.0=P.0, control=list(fnscale=-1, trace=trace));
  if (trace) cat("Finished optim.\n")

  ## optim.out = optim(abphi.0, blogit.llh.mm, gr=grad.blogit.llh.mm, method="BFGS", hessian=TRUE,
  ##   y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate, m.0=m.0, P.0=P.0, control=list(fnscale=-1));
  
  ## blogit.llh.mm(abphi.0, y, X.re, X.fe, n, shape, rate, m.0, P.0)
  ## grad.blogit.llh.mm(abphi.0, y, X.re, X.fe, n, shape, rate, m.0, P.0)
  ## hessian.blogit.llh.mm(abphi.0, y, X.re, X.fe, n, shape, rate, P.0)

  abphi.pm = optim.out$par;

  ## num - numerically
  ## pen - analytically (pen and paper)
  
  ## hess = optim.out$hessian
  hess.pen = hessian.blogit.llh.mm(abphi.pm, y, X.re, X.fe, n, shape, rate, P.0)
  hess = hess.pen
  
  ## V.num = solve(-1*optim.out$hessian); ## I encountered some numerical instability in the German dataset.
  ## V.pen = solve(-1 * hess.pen);
  ## Chol.ml = chol(V.pen);
  ## V = V.pen;
  
  ## m = as.numeric(beta.pm);

  V = NA;
  if (calc.V) V = solve(-1 * hess);

  out = list("V"=V, "H"=hess, "H.pen"=hess.pen, "convergence"=optim.out$convergence, "m"=abphi.pm)
}

blogit.mm.map.alt <- function(y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                              m.0=array(0, dim=c(ncol(X.fe))),
                              P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                              abphim.0=NULL, calc.V=TRUE, trace=FALSE)
{ 
  X = cbind(X.re, X.fe);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(abphim.0)) { abphim.0 = matrix(0, P+2); abphim.0[P+1] = 1; }

  optim.out <- optim(abphim.0, blogit.llh.mm.alt, hessian=TRUE,
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate, kappa=kappa,
                     m.0=m.0, P.0=P.0, control=list(fnscale=-1, trace=trace));

  abphim.pm = optim.out$par;
  hess = optim.out$hessian

  V = NA;
  if (calc.V) V = solve(-1 * hess);

  out = list("V"=V, "H"=hess, "convergence"=optim.out$convergence, "m"=abphim.pm)
}

################################################################################

################################################################################

llh.norm <- function(x, kappa, UorP, is.prec=TRUE, is.mean=TRUE)
{
  ## x: point at which llh is evaluated.
  ## kappa: Prec * mean;
  ## UorP: chol(Prec) or Prec
  ## is.prec: specifies that UorP is Prec.
  
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

################################################################################
                              ## IND. METROP 1 ##
################################################################################

log.tget.to.ppsl.blogit <- function(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0, m, UorP, df=Inf, is.prec=TRUE)
{
  ## Calculate log f - log q
  if (df==Inf)
    log.fdivq = llh.norm(abphi, m, UorP, is.prec)
  else
    log.fdivq = llh.t(abphi, m, UorP, df, is.prec)

  log.fdivq = blogit.llh.mm(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0) - log.fdivq
  
  log.fdivq
}

draw.abphi.ind.MH <- function(abphi, y, X.re, X.fe, n, shape, rate, m0, P0, map, U.map, log.fdivq, df=Inf)
{
  ## Ind MH sample when m0, P0 might be changing.
  ## beta: previous beta
  ## y, X, n: data
  ## m0, P0: priors for beta
  ## map: posterior mode
  ## U.map: chol(-H) where H is Hessian at mode
  ## df: degrees of freedom for proposal
  ## log.fdiq.q : previous

  X = cbind(X.re, X.fe)
  
  P.a  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.ab = ncol(X);
  
  P = length(abphi)
  zero = rep(0, P);

  ## log.fdivq = log.tget.to.ppsl(beta, y, X, n, m.0, P.0, map, U.map, df=df, is.prec=FALSE)

  ## Propose
  if (df==Inf) ep.ppsl = rnorm(P) else ep.ppsl = rt(P, df=df)

  ppsl = map + backsolve(U.map, ep.ppsl);

  log.fdivq.ppsl = log.tget.to.ppsl.blogit(ppsl, y, X.re, X.fe, n, shape, rate, m0, P0, map, U.map, df=df, is.prec=FALSE)

  ## acceptance prob
  a.prob = min( exp(log.fdivq.ppsl - log.fdivq), 1);
  accept = runif(1) < a.prob

  if (accept) {
    abphi = ppsl;
    log.fdivq = log.fdivq.ppsl
  }

  out = list("abphi"=abphi, "a.prob"=a.prob, "accept"=accept, "log.fdivq"=log.fdivq)

  out
}

ind.metropolis.blogit <- function(y, X.re, X.fe, n=rep(1, length(y)), shape=1, rate=1,
                                  m.0=rep(0, ncol(X.fe)), P.0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe)),
                                  samp=1000, burn=100, verbose=1000, df=Inf)
{
  ## An independent MH routine that can be adapted to varying P0.
  ## Assume proper posterior
  
  ## y: response, number of responses in each category
  ## X: design
  ## n: number of trials per draw.
  
  y = as.matrix(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  
  N = nrow(X.re);
  P.a  = ncol(X.re); 
  P.b  = ncol(X.fe); 
  P.ab = P.a + P.b

  a.idc = 1:P.a;
  b.idc = 1:P.b + P.a;
  phi.idc = P.ab + 1

  ## n = rep(1, N)
  
  if (is.null(m.0)) m.0 = rep(0, P.b);
  if (is.null(P.0)) {
    p.0 = 0e-4
    P.0 = diag(p.0, P.b);
  }

  out <- list(beta = array(0, dim=c(samp, P.b)),
              alpha = array(0, dim=c(samp, P.a)),
              ab = array(0, dim=c(samp, P.ab)),
              phi = rep(0, samp),
              a.prob = rep(0, samp)
              )

  ## Find mode
  out.map = blogit.mm.map(y, X.re, X.fe, n, shape, rate, m.0, P.0, abphi.0=NULL)
  m = out.map$m
  U = chol(-1 * out.map$H)
  cat("Finished with MLE.  Convergence =", out.map$convergence, "\n")

  ## start at ppsl mean.
  ep   = rep(0, P.ab+1)
  abphi = out.map$m
  log.fdivq = log.tget.to.ppsl.blogit(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0, m, U, df=df, is.prec=FALSE)
  
  ## Generate proposal
  out.abphi = draw.abphi.ind.MH(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0, map=m, U.map=U, log.fdivq=log.fdivq, df=df)
  abphi = out.abphi$abphi
  log.fdivq = out.abphi$log.fdivq

  ## Timing
  start.time = proc.time()
  naccept = 0

  ## Do Metropolis ##
  for (i in 1:(samp+burn)) {

    if (i==burn+1) { start.ess = proc.time(); }

    ## Generate proposal
    out.abphi = draw.abphi.ind.MH(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0, map=m, U.map=U, log.fdivq=log.fdivq, df=df)
    abphi = out.abphi$abphi
    log.fdivq = out.abphi$log.fdivq
    
    if (i > burn) {
      out$beta[i-burn,]  = out.abphi$abphi[b.idc]
      out$alpha[i-burn,] = out.abphi$abphi[a.idc]
      out$ab[i-burn,]    = out.abphi$abphi[1:P.ab]
      out$phi[i-burn]    = out.abphi$abphi[phi.idc]
      out$a.prob[i-burn] = out.abphi$a.prob
      naccept = naccept + out.abphi$accept
    }

    if (i %% verbose == 0) {
      if (i > burn) cat("Ind MH MM: Ave a.prob:", mean(out$a.prob[1:(i-burn)]), ", ");
      cat("Iteration:", i, "\n");
    }
    
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$naccept    = naccept
  out$acceptr    = naccept / samp;
  out$m          = out.map$m
  out$V          = out.map$V

  out
}

################################################################################
                              ## IND. METROP 2 ##
################################################################################

log.tget.to.ppsl.blogit.2 <- function(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0,
                                      m.ab, U.ab, shape.ppsl, rate.ppsl, df=Inf)
{
  P.ab = length(abphi) - 1;
  
  ab  = abphi[1:P.ab]
  phi = abphi[P.ab+1]
  
  ## Calculate log f - log q
  if (df==Inf)
    log.fdivq = llh.norm(ab, m.ab, U.ab, FALSE)
  else
    log.fdivq = llh.t(ab, m.ab, U.ab, df, FALSE)

  log.fdivq = log.fdivq + dgamma(phi, shape=shape.ppsl, rate=rate.ppsl, log=TRUE);

  log.fdivq = blogit.llh.mm(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0) - log.fdivq
  
  log.fdivq
}

draw.abphi.ind.MH.2 <- function(abphi, y, X.re, X.fe, n, shape, rate, m0, P0,
                                m.ab, U.ab, shape.ppsl, rate.ppsl, log.fdivq, df=Inf)
{
  ## Ind MH sample when m0, P0 might be changing.
  ## beta: previous beta
  ## y, X, n: data
  ## m0, P0: priors for beta
  ## map: posterior mode
  ## U.map: chol(-H) where H is Hessian at mode
  ## df: degrees of freedom for proposal
  ## log.fdiq.q : previous

  X = cbind(X.re, X.fe)
  
  P.a  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.ab = ncol(X);
  
  P = length(abphi)
  zero = rep(0, P);

  ## log.fdivq = log.tget.to.ppsl(beta, y, X, n, m.0, P.0, map, U.map, df=df, is.prec=FALSE)

  ## Propose
  if (df==Inf) ep.ppsl = rnorm(P.ab) else ep.ppsl = rt(P, df=df)

  ab.ppsl  = m.ab + backsolve(U.ab, ep.ppsl);
  phi.ppsl = rgamma(1, shape.ppsl, rate=rate.ppsl)

  ppsl = c(ab.ppsl, phi.ppsl)

  log.fdivq.ppsl = log.tget.to.ppsl.blogit.2(ppsl, y, X.re, X.fe, n, shape, rate, m0, P0, m.ab, U.ab, shape.ppsl, rate.ppsl, df=df)

  ## acceptance prob
  a.prob = min( exp(log.fdivq.ppsl - log.fdivq), 1);
  accept = runif(1) < a.prob

  if (accept) {
    abphi = ppsl;
    log.fdivq = log.fdivq.ppsl
  }

  out = list("abphi"=abphi, "a.prob"=a.prob, "accept"=accept, "log.fdivq"=log.fdivq)

  out
}

ind.metropolis.blogit.2 <- function(y, X.re, X.fe, n=rep(1, length(y)), shape=1, rate=1,
                                    m.0=rep(0, ncol(X.fe)), P.0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe)),
                                    samp=1000, burn=100, verbose=1000, df=Inf)
{
  ## An independent MH routine that can be adapted to varying P0.
  ## Assume proper posterior
  
  ## y: response, number of responses in each category
  ## X: design
  ## n: number of trials per draw.
  
  y = as.matrix(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  
  N = nrow(X.re);
  P.a  = ncol(X.re); 
  P.b  = ncol(X.fe); 
  P.ab = P.a + P.b

  a.idc = 1:P.a;
  b.idc = 1:P.b + P.a;
  phi.idc = P.ab + 1
  ab.idc = 1:P.ab

  ## n = rep(1, N)
  
  if (is.null(m.0)) m.0 = rep(0, P.b);
  if (is.null(P.0)) {
    p.0 = 0e-4
    P.0 = diag(p.0, P.b);
  }

  out <- list(beta = array(0, dim=c(samp, P.b)),
              alpha = array(0, dim=c(samp, P.a)),
              ab = array(0, dim=c(samp, P.ab)),
              phi = rep(0, samp),
              a.prob = rep(0, samp)
              )

  ## Find mode
  out.map = blogit.mm.map(y, X.re, X.fe, n, shape, rate, m.0, P.0, abphi.0=NULL)
  m = out.map$m
  cat("Finished with MLE.  Convergence =", out.map$convergence, "\n")

  ## Setup necessary quantities, mean, chol(Prec), shape.ppsl, rate.ppsl
  m.ab = m[ab.idc];
  V.ab = out.map$V[ab.idc, ab.idc]
  U.ab = chol( solve(V.ab) );

  m.phi = m[phi.idc]
  V.phi = out.map$V[phi.idc, phi.idc]
  shape.ppsl = m.phi^2 / V.phi
  rate.ppsl  = m.phi   / V.phi
  
  ## start at ppsl mean.
  ep   = rep(0, P.ab+1)
  abphi = out.map$m
  log.fdivq <- log.tget.to.ppsl.blogit.2(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0,
                                         m.ab, U.ab, shape.ppsl, rate.ppsl, df=df)
  
  ## Generate proposal
  out.abphi <- draw.abphi.ind.MH.2(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0,
                                   m.ab, U.ab, shape.ppsl, rate.ppsl, log.fdivq=log.fdivq, df=df)
  abphi = out.abphi$abphi
  log.fdivq = out.abphi$log.fdivq

  ## Timing
  start.time = proc.time()
  naccept = 0

  ## Do Metropolis ##
  for (i in 1:(samp+burn)) {

    if (i==burn+1) { start.ess = proc.time(); }

    ## Generate proposal
    out.abphi <- draw.abphi.ind.MH.2(abphi, y, X.re, X.fe, n, shape, rate, m.0, P.0,
                                     m.ab, U.ab, shape.ppsl, rate.ppsl, log.fdivq=log.fdivq, df=df)
    abphi = out.abphi$abphi
    log.fdivq = out.abphi$log.fdivq
    
    if (i > burn) {
      out$beta[i-burn,]  = out.abphi$abphi[b.idc]
      out$alpha[i-burn,] = out.abphi$abphi[a.idc]
      out$ab[i-burn,]    = out.abphi$abphi[1:P.ab]
      out$phi[i-burn]    = out.abphi$abphi[phi.idc]
      ## out$abphi[i-burn,] = abphi
      out$a.prob[i-burn] = out.abphi$a.prob
      naccept = naccept + out.abphi$accept
    }

    if (i %% verbose == 0) {
      if (i > burn) cat("Ind MH MM 2: Ave a.prob:", mean(out$a.prob[1:(i-burn)]), ", ");
      cat("Iteration:", i, "\n");
    }
    
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$naccept    = naccept
  out$acceptr    = naccept / samp;
  out$m.ab       = m.ab
  out$V.ab       = V.ab
  out$m.phi      = m.phi
  out$V.phi      = V.phi

  out
}

################################################################################
                                    ## PG ##
################################################################################

## Bayesian logistic regression
##------------------------------------------------------------------------------
logit.PG.mm <- function(y, X.re, X.fe, n=rep(1, length(y)), shape=1, rate=1,
                        m0=rep(0, ncol(X.fe)), P0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe)),
                        samp=1000, burn=500, verbose=500)
{
  ## X: n by p matrix
  ## y: n by 1 vector, total # response
  ## n: n by 1 vector, # of obs at distinct x

  y = as.numeric(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  X    = cbind(X.re, X.fe)
  
  P.a  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.ab = ncol(X)
  N = nrow(X)
  
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a

  Z = colSums(X * (y-n/2));
  Z[b.idc] = Z[b.idc] + P0 %*% m0;

  ## PsiToBeta = solve(t(X) %*% X) %*% t(X);

  w = rep(0,N)
  ## w = w.known;
  ab = rep(0.0, P.ab)
  phi = shape / rate;

  output <- list(w = matrix(nrow=samp, ncol=N),
                 ab = matrix(nrow=samp, ncol=P.ab),
                 alpha = matrix(nrow=samp, ncol=P.a),
                 beta  = matrix(nrow=samp, ncol=P.b),
                 phi = rep(0, samp)
                 )

  ## c_k = (1:200-1/2)^2 * pi^2 * 4;

  ## Timing
  start.time = proc.time()
  
  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();
    
    ## draw w
    psi = drop(X%*%ab)
    
    ## Devroye is faster anyway.
    w = rpg.devroye(N, n, psi);
    
    # draw beta - Joint Sample.
    PP = t(X) %*% (X * w);
    PP[b.idc, b.idc] = PP[b.idc,b.idc] + P0
    PP[a.idc, a.idc] = PP[a.idc, a.idc] + diag(phi, P.a);
    ## U = chol(PP);
    ## m = backsolve(U, Z, transpose=TRUE);
    ## m = backsolve(U, m);
    ## ab = m + backsolve(U, rnorm(P.ab))
    S = chol2inv(chol(PP));
    m = S %*% as.vector(Z);
    ab = m + t(chol(S)) %*% rnorm(P.ab);
    alpha= ab[a.idc]
    
    phi = rgamma(1, shape=0.5*(shape+P.a), rate=0.5*(t(alpha) %*% alpha + rate));
    
    # Record if we are past burn-in.
    if (j>burn) {
        output$w[j-burn,] <- w
        output$ab[j-burn,] <- ab
        output$alpha[j-burn,] <- alpha
        output$beta[j-burn,] <- ab[b.idc]
        output$phi[j-burn] <- phi
    }

    if (j %% verbose == 0) { print(paste("LogitPG MM: Iteration", j)); }
  }

  end.time = proc.time()
  output$total.time = end.time - start.time
  output$ess.time   = end.time - start.ess

  ## ## Add new data to output.
  ## output$"y" = y;
  ## output$"X" = X;
  ## output$"n" = n;

  output
} ## logit.gibbs.R

################################################################################
                                   ## TEST ##
################################################################################

if (FALSE) {

  N = 500;
  P = 1;
  groups = 5;

  group.id = matrix(1:groups, ncol=groups, nrow=ceiling(N / groups), byrow=TRUE)[1:N];
  group.id = factor(group.id);
  X.re = as.matrix( model.matrix(~ group.id + 0) )
  
  ## Correlated predictors
  rho = 0.5
  Sig = matrix(rho, nrow=P, ncol=P); diag(Sig) = 1.0;
  U   = chol(Sig);
  X.fe = matrix(rnorm(N*P), nrow=N, ncol=P) %*% U;

  X = cbind(X.re, X.fe);

  phi = 1

  alpha = rnorm(groups, 0, sqrt(1/phi));
  beta = rnorm(P, mean=0, sd=2);
  
  psi = X %*% c(alpha, beta);
  p = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);
  n = rep(1, N);

}

if (FALSE) {
  
  P.a = ncol(X.re);
  P.b = ncol(X.fe);

  shape = 1
  rate  = 1
  
  prec.a = rep((shape-1) / rate, P.a)
  prec.b = rep(0.01, P.b)
  
  m.0   = rep(0, P.b);
  P.0   = diag(prec.b, P.b);
  
  m.0.mlogit = matrix(0, nrow=P.a+P.b, ncol=1);
  P.0.mlogit = array(diag(c(prec.a, prec.b), P.a+P.b), dim=c(P.a+P.b, P.a+P.b, 1))

  m.0.PG.mm = m.0.mlogit[,1];
  P.0.PG.mm = P.0.mlogit[,,1];
  
  out.map = blogit.mm.map(y, X.re, X.fe, n, shape, rate, m.0, P.0, abphi.0=NULL)
  m = out.map$m
  U = chol(-1 * out.map$H)
  V = chol2inv(U)

  samp = 10000
  burn = 1000
  verbose = 1000
  df = Inf
  out.ind = ind.metropolis.blogit(y, X.re, X.fe, n, shape, rate, m.0, P.0, samp=samp, burn=burn, verbose=verbose, df=df)
  out.ind.2 = ind.metropolis.blogit.2(y, X.re, X.fe, n, shape, rate, m.0, P.0, samp=samp, burn=burn, verbose=verbose, df=df)

  out.pg.mm = logit.PG.mm(y, X.re, X.fe, n, shape, rate, m.0, P.0, samp=samp, burn=burn, verbose=verbose);

  out.mlogit = mlogit.MH.R(y, X, n, m.0=m.0.mlogit, P.0=P.0.mlogit, beta.0=NULL,
    samp=samp, burn=burn, method="Ind", tune=1.0, df=df, verbose=1000)
  
  out.pg = logit.R(y, X, m0=m.0.PG.mm, P0=P.0.PG.mm, samp=samp, burn=burn, verbose=verbose);
  
  rbind(sum.stat(out.ind$ab, out.ind$ess.time[3]),
        sum.stat(out.ind$phi, out.ind$ess.time[3]))

  rbind(sum.stat(out.ind.2$ab, out.ind$ess.time[3]),
        sum.stat(out.ind.2$phi, out.ind$ess.time[3]))

  rbind(rbind(sum.stat(out.pg.mm$ab, out.pg.mm$ess.time[3])),
        rbind(sum.stat(out.pg.mm$phi, out.pg.mm$ess.time[3])))

  rbind(sum.stat(out.pg$beta, out.pg$ess.time[3]))

  rbind(sum.stat(out.mlogit$beta, out.mlogit$ess.time[3]))

  sstat.ind.1 = sum.stat(out.ind$ab, out.ind$ess.time[3]);
  sstat.ind.2 = sum.stat(out.ind.2$ab, out.ind$ess.time[3]);
  sstat.pg.mm = sum.stat(out.pg.mm$ab, out.pg.mm$ess.time[3]);
  
}
