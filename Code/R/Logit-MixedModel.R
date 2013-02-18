
## binary logit with groups with different intercepts.

## psi_i = m + delta_j + x_i beta; alpha_j = m + delta_j; i in group j.

## 1: alpha, beta, phi; normal proposal
## 2: alpha, beta, phi; normal ppsl for alpha, beta, gamma ppsl for phi
## 3: alpha, beta, m, theta; theta = log(phi); normal ppsl
## 4: alpha, beta, m, phi; normal ppsl
## 6: delta, beta, m, theta; theta = log(phi); normal ppsl
## 7: delta, beta, m, phi; normal ppsl

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
  hess[b.idc, b.idc] = hess[b.idc, b.idc] - P.0
  hess[a.idc, a.idc] = hess[a.idc, a.idc] - diag(phi, P.a);
  hess[a.idc, phi.idc] = -1 * alpha;
  hess[phi.idc, a.idc] = -1 * alpha;
  hess[phi.idc, phi.idc] = -1 * (0.5 * (P.a + shape) - 1) / phi^2;

  hess
}

blogit.llh.mm.3 <- function(abtm, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))))
{
  ## abtm: alpha, beta, theta, m
  ## \psi_t = X.fe alpha + X.re \beta
  ## \phi = log \theta
  ## alpha \sim N(m, 1/\phi)
  ## m \sim N(0, \kappa^2 / \phi)
  
  ## abphi: (P.a + P.b + 1) x J-1; coefficients, column assumed to be beta_J = 0.
  ## y : N; number of responses
  ## X : N x P.ab: design matrix
  ## P.0 : P.b x P.b; array of matrices for independent prior. 
  
  N    = length(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  X    = cbind(X.re, X.fe)
  P.a  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.ab = P.a + P.b
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  t.idc = P.ab + 1;

  alpha = abtm[a.idc]
  beta  = abtm[b.idc]
  theta = abtm[t.idc]
  if (kappa != 0) m = abtm[t.idc+1] else m = 0;
  not.flat = as.numeric(kappa != Inf);
  
  ab = abtm[1:P.ab]
  Psi = X %*% ab;
  ## p = exp(Psi) / (1 + exp(Psi));

  ## Based on p(y|alpha, beta) p(beta) p(alpha|theta) p(theta)
  llh = t(y) %*% Psi - sum(n * log(1 + exp(Psi)));
  llh = llh - 0.5 * t(beta - m.0) %*% P.0 %*% (beta - m.0); 
  llh = llh + 0.5 * (P.a + shape) * theta - 0.5 * (sum((alpha-m)^2) + rate) * exp(theta); 
  if (kappa != 0) llh = llh + (0.5 * theta) * not.flat - 0.5 * m^2 / kappa^2 * exp(theta); 

  llh
}

hessian.blogit.mm.3 <- function(abtm, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                                P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))) )

{
  ## Different parameterization of prior, log gamma, produces different posteriors.
  N    = length(y)
  X    = cbind(X.re, X.fe);
  P.a  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.ab = P.a + P.b
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  t.idc = P.ab + 1;
  m.idc = t.idc + 1;
  P.all = P.a + P.b + 1 + as.numeric(kappa!=0)

  alpha = abtm[a.idc]
  beta  = abtm[b.idc]
  theta = abtm[t.idc]
  
  ab = abtm[1:P.ab]
  Psi = X %*% ab;
  p   = exp(Psi) / (1 + exp(Psi));
  phi = exp(theta);

  if (kappa != 0) m = abtm[m.idc] else m = 0;

  d = (p^2 - p) * n
  H = t(X) %*% (X * rep(d, P.ab));
  
  hess = matrix(0, P.all, P.all);
  hess[1:P.ab, 1:P.ab] = H
  hess[b.idc, b.idc] = hess[b.idc, b.idc] - P.0
  hess[a.idc, a.idc] = hess[a.idc, a.idc] - diag(phi, P.a)
  hess[a.idc, t.idc] = -1 * (alpha-m) * phi
  hess[t.idc, a.idc] = -1 * (alpha-m) * phi
  hess[t.idc, t.idc] = -0.5 * (t(alpha-m) %*% (alpha-m) + rate) * phi
  if (kappa != 0) {
    hess[m.idc, m.idc] = -1 / kappa^2 * phi - P.a * phi
    hess[m.idc, t.idc] = -m / kappa^2 * phi + sum(alpha-m) * phi
    hess[t.idc, m.idc] = hess[m.idc, t.idc]
    hess[m.idc, a.idc] = phi
    hess[a.idc, m.idc] = phi
    hess[t.idc, t.idc] = hess[t.idc, t.idc] - 0.5 * m^2 / kappa^2 * phi
  }

  hess
}

blogit.llh.mm.4 <- function(abphim, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))))                         
{
  ## abphim: alpha, beta, phi, m
  ## \psi_t = X.fe alpha + X.re \beta
  ## \phi = log \theta
  ## alpha \sim N(m, 1/\phi)
  ## m \sim N(0, \kappa^2 / \phi)
  
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
  beta  = abphim[b.idc]
  phi   = abphim[phi.idc]
  if (kappa != 0) m = abphim[phi.idc+1] else m = 0;
  not.flat = as.numeric(kappa != Inf);

  if (phi < 0) return(-Inf)
  
  ab = abphim[1:P.ab]
  Psi = X %*% ab;
  ## p = exp(Psi) / (1 + exp(Psi));

  llh = t(y) %*% Psi - sum(n * log(1 + exp(Psi)));
  llh = llh - 0.5 * t(beta - m.0) %*% P.0 %*% (beta - m.0);
  llh = llh + (0.5 * (P.a + shape) - 1) * log(phi) - 0.5 * phi * ( sum((alpha-m)^2) + rate );
  if (kappa != 0) llh = llh + 0.5 * log(phi) * not.flat - 0.5 * phi * m^2/kappa^2;

  llh
}

## Need to check
grad.blogit.mm.4 <- function(abphim, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
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
  f.idc = P.ab + 1;
  m.idc = f.idc + 1;

  alpha = abphim[a.idc]
  beta  = abphim[b.idc]
  phi   = abphim[f.idc]
  inc.m = as.numeric(kappa!=0)
  if (inc.m) m = abphim[f.idc+1] else m = 0;
  not.flat = as.numeric(kappa != Inf);
  
  if (phi < 0) return(NA)
  
  ab = abphim[1:P.ab]
  Psi = X %*% ab;
  p = exp(Psi) / (1 + exp(Psi));

  grad = rep(0, P.ab+1+inc.m)
  grad[1:P.ab] = t(y - p * n) %*% X;
  grad[b.idc]  = grad[b.idc] - P.0 %*% (beta - m.0)
  grad[a.idc]  = grad[a.idc] - phi * (alpha - m);
  grad[f.idc]  = (0.5 * (P.a + shape) - 1 ) / phi - 0.5 * (t(alpha-m) %*% (alpha-m) + rate);
  if (kappa!=0) {
    grad[f.idc] = grad[f.idc] + (0.5 / phi) * not.flat - 0.5 * m^2 / kappa^2;
    grad[m.idc] = (sum(alpha-m) - m / kappa^2) * phi
  }
  
  grad
}

hessian.blogit.mm.4 <- function(abphim, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                                P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))) )
{  
  ## Using phi parameterization
  N    = length(y)
  X    = cbind(X.re, X.fe);
  P.a  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.ab = P.a + P.b
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  f.idc = P.ab + 1;
  m.idc = f.idc + 1;
  P.all = P.a + P.b + 1 + as.numeric(kappa!=0)

  alpha = abphim[a.idc]
  beta  = abphim[b.idc]
  phi   = abphim[f.idc]
  
  ab = abphim[1:P.ab]
  Psi = X %*% ab;
  p   = exp(Psi) / (1 + exp(Psi));

  if (kappa != 0) m = abphim[m.idc] else m = 0.0;
  not.flat = as.numeric(kappa != Inf);

  d = (p^2 - p) * n
  H = t(X) %*% (X * rep(d, P.ab));
  
  hess = matrix(0, P.all, P.all);
  hess[1:P.ab, 1:P.ab] = H
  hess[b.idc, b.idc] = hess[b.idc, b.idc] - P.0
  hess[a.idc, a.idc] = hess[a.idc, a.idc] - diag(phi, P.a)
  hess[a.idc, f.idc] = -1 * (alpha-m) 
  hess[f.idc, a.idc] = -1 * (alpha-m) 
  hess[f.idc, f.idc] = -1 * (0.5 * (P.a + shape) - 1 ) / phi^2
  if (kappa != 0) {
    hess[m.idc, m.idc] = -1 / kappa^2 * phi - P.a * phi
    hess[m.idc, f.idc] = -m / kappa^2 + sum(alpha-m)
    hess[f.idc, m.idc] = hess[m.idc, f.idc]
    hess[m.idc, a.idc] = phi
    hess[a.idc, m.idc] = phi
    hess[f.idc, f.idc] = hess[f.idc, f.idc] - (0.5 / phi^2) * not.flat
  }

  hess
}

blogit.llh.mm.5 <- function(abtm, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))))                         
{
  ## THIS IS INCORRECT.  USE VERSION 3.  REPARAMETERIZATION SHOULD CHANGE POSTERIOR.
  
  ## abtm: alpha, beta, theta, m
  ## \psi_t = X.fe alpha + X.re \beta
  ## \phi = log \theta
  ## alpha \sim N(m, 1/\phi)
  ## m \sim N(0, \kappa^2 / \phi)
  
  ## this finds the mode under a different parameterization.  But this should
  ## not be used for sampling phi, the precision.

  ## N    = length(y)
  ## X.re = as.matrix(X.re)
  ## X.fe = as.matrix(X.fe)
  ## X    = cbind(X.re, X.fe)
  ## P.a  = ncol(X.re)
  ## P.b  = ncol(X.fe)
  ## P.ab = P.a + P.b
  ## a.idc = 1:P.a
  ## b.idc = 1:P.b + P.a
  ## t.idc = P.ab + 1;

  ## alpha = abtm[a.idc]
  ## beta  = abtm[b.idc]
  ## theta = abtm[t.idc]
  ## if (kappa != 0) m = abtm[t.idc+1] else m = 0;
  
  ## ab = abtm[1:P.ab]
  ## Psi = X %*% ab;
  ## ## p = exp(Psi) / (1 + exp(Psi));

  ## ## Based on p(y|alpha, beta) p(beta) p(alpha|theta) p(theta)
  ## llh = t(y) %*% Psi - sum(n * log(1 + exp(Psi)));
  ## llh = llh - 0.5 * t(beta - m.0) %*% P.0 %*% (beta - m.0); 
  ## llh = llh + (0.5 * (P.a + shape) - 1) * theta - 0.5 * (sum((alpha-m)^2) + rate) * exp(theta); 
  ## if (kappa != 0) llh = llh + 0.5 * theta - 0.5 * m^2 / kappa^2 * exp(theta); 

  llh = blogit.llh.mm.3(abtm, y, X.re, X.fe, n, shape=shape-2, rate, kappa, m.0, P.0);
  
  llh
}

hessian.blogit.mm.5 <- function(abtm, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                                  ## m.0=array(0, dim=c(ncol(X.fe))),
                                  P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))) )

{
  ## log phi = theta so that term drops out of hessian.  Hence hessian.3 = hessian.5
  hess = hessian.blogit.mm.3(abtm, y, X.re, X.fe, n, shape, rate, kappa, P.0)
}

blogit.llh.mm.6 <- function(dbmt, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))))
{
  ## dbmp: delta, beta, m, phi
  ## \psi_t = X.fe (m + \delta) + X.re \beta
  ## \phi = log \theta
  ## delta \sim N(0, 1/\phi)
  ## m \sim N(0, \kappa^2 / \phi)
  if (kappa == 0) return(NA);
  
  N    = length(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  X    = cbind(X.re, X.fe, 1)
  P.d  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.dbm = P.d + P.b + 1
  d.idc = 1:P.d
  b.idc = 1:P.b + P.d
  m.idc = P.dbm
  t.idc = P.dbm + 1;
  
  delta = dbmt[d.idc]
  beta  = dbmt[b.idc]
  m     = dbmt[m.idc]
  theta = dbmt[t.idc]
  not.flat = as.numeric(kappa != Inf);

  ## print(dbmt)
  ## print(delta)
  
  dbm = dbmt[1:P.dbm]
  Psi = X %*% dbm
  ## p = exp(Psi) / (1 + exp(Psi));

  ## Based on p(y|delta, beta) p(beta) p(delta|theta) p(theta)
  llh = t(y) %*% Psi - sum(n * log(1 + exp(Psi)));
  llh = llh - 0.5 * t(beta - m.0) %*% P.0 %*% (beta - m.0); 
  llh = llh + 0.5 * (P.d + shape) * theta - 0.5 * (t(delta) %*% delta + rate) * exp(theta); 
  llh = llh + (0.5 * theta) * not.flat - 0.5 * m^2 / kappa^2 * exp(theta); 

  llh
}

hessian.blogit.mm.6 <- function(dbmt, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                                P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))) )
{
  ## Different parameterization of prior, log gamma, produces different posteriors.
  N    = length(y)
  X    = cbind(X.re, X.fe, 1);
  P.d  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.dbm = P.d + P.b + 1
  d.idc = 1:P.d
  b.idc = 1:P.b + P.d
  m.idc = P.dbm
  t.idc = P.dbm + 1;
  P.all = P.dbm + 1
 
  delta = dbmt[d.idc]
  m     = dbmt[m.idc]
  beta  = dbmt[b.idc]
  theta = dbmt[t.idc]

  dbm = dbmt[1:P.dbm]
  Psi = cbind(X) %*% dbm;
  p   = exp(Psi) / (1 + exp(Psi));
  phi = exp(theta);

  dg = (p^2 - p) * n
  H = t(X) %*% (X * rep(dg, P.dbm));
  
  hess = matrix(0, P.all, P.all);
  hess[1:P.dbm, 1:P.dbm] = H
  hess[b.idc, b.idc] = hess[b.idc, b.idc] - P.0
  hess[d.idc, d.idc] = hess[d.idc, d.idc] - diag(phi, P.d)
  hess[d.idc, t.idc] = -1 * (delta) * phi
  hess[t.idc, d.idc] = -1 * (delta) * phi
  hess[t.idc, t.idc] = -0.5 * (t(delta) %*% (delta) + rate) * phi - 0.5 * m^2 / kappa^2 * phi
  hess[m.idc, m.idc] = hess[m.idc, m.idc] -1 / kappa^2 * phi
  hess[m.idc, t.idc] = -m / kappa^2 * phi
  hess[t.idc, m.idc] = hess[m.idc, t.idc]

  hess
}

blogit.llh.mm.7 <- function(dbmp, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))))                         
{
  ## dbmp: delta, beta, m, phi
  ## \psi_t = X.fe (m + \delta) + X.re \beta
  ## delta \sim N(0, 1/\phi)
  ## m \sim N(0, \kappa^2 / \phi)
  if (kappa == 0) return(NA);
  
  N = length(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  X   = cbind(X.re, X.fe, 1)
  P.d = ncol(X.re)
  P.b = ncol(X.fe)
  P.dbm = P.d + P.b + 1
  d.idc = 1:P.d
  b.idc = 1:P.b + P.d
  m.idc = P.dbm
  p.idc = P.dbm + 1;

  delta = dbmp[d.idc]
  m     = dbmp[m.idc]
  beta  = dbmp[b.idc]
  phi   = dbmp[p.idc]
  not.flat = as.numeric(kappa != Inf);

  if (phi < 0) return(-Inf)
  
  dbm = dbmp[1:P.dbm]
  Psi = X %*% dbm;
  ## p = exp(Psi) / (1 + exp(Psi));

  llh = t(y) %*% Psi - sum(n * log(1 + exp(Psi)));
  llh = llh - 0.5 * t(beta - m.0) %*% P.0 %*% (beta - m.0);
  llh = llh + (0.5 * (P.d + shape) - 1) * log(phi) - 0.5 * phi * ( sum((delta)^2) + rate );
  llh = llh + 0.5 * log(phi) * not.flat - 0.5 * phi * m^2/kappa^2;

  llh
}

hessian.blogit.mm.7 <- function(dbmp, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                                P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))) )
{
  ## Using phi parameterization
  N    = length(y)
  X    = cbind(X.re, X.fe, 1);
  P.d  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.dbm = P.d + P.b + 1
  d.idc = 1:P.d
  b.idc = 1:P.b + P.d
  m.idc = P.dbm;
  p.idc = P.dbm + 1;
  P.all = P.dbm + 1;
    
  delta = dbmp[d.idc]
  m     = dbmp[m.idc]
  beta  = dbmp[b.idc]
  phi   = dbmp[p.idc]
  
  dbm = dbmp[1:P.dbm]
  Psi = X %*% dbm;
  p   = exp(Psi) / (1 + exp(Psi));
  not.flat = as.numeric(kappa != Inf);

  dg = (p^2 - p) * n
  H = t(X) %*% (X * rep(dg, P.dbm));
  
  hess = matrix(0, P.all, P.all);
  hess[1:P.dbm, 1:P.dbm] = H
  hess[b.idc, b.idc] = hess[b.idc, b.idc] - P.0
  hess[d.idc, d.idc] = hess[d.idc, d.idc] - diag(phi, P.d)
  hess[d.idc, p.idc] = -1 * (delta) 
  hess[p.idc, d.idc] = -1 * (delta) 
  hess[p.idc, p.idc] = -1 * (0.5 * (P.d + shape) - 1 ) / phi^2 - (0.5 / phi^2) * not.flat
  hess[m.idc, m.idc] = hess[m.idc, m.idc] + -1 / kappa^2 * phi
  hess[m.idc, p.idc] = -m / kappa^2 
  hess[p.idc, m.idc] = hess[m.idc, p.idc]

  hess
}

blogit.llh.mm.8 <- function(ab, y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))))                         
{
  ## abphim: alpha, beta, phi, m
  ## \psi_t = X.fe alpha + X.re \beta
  ## \phi = log \theta
  ## alpha \sim N(m, 1/\phi)
  ## m \sim N(0, \kappa^2 / \phi)
  ## Now integrate out phi and m.
  ## \alpha sim N(0, (1+\kappa^2) / \phi)
  ## \alpha sim t(shape, 0, rate * (1+\kappa)^2).
  ## I need to check the last equation.
  
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

  alpha = ab[a.idc]
  beta  = ab[b.idc]

  ab = ab[1:P.ab]
  Psi = X %*% ab;
  ## p = exp(Psi) / (1 + exp(Psi));

  llh = t(y) %*% Psi - sum(n * log(1 + exp(Psi)));
  llh = llh - 0.5 * t(beta - m.0) %*% P.0 %*% (beta - m.0);
  llh = llh - 0.5 * (shape + P.a) * log( 1 + t(alpha) %*% alpha / (rate * (1 + kappa^2)) )
  llh
}

################################################################################
##------------------------------------------------------------------------------

blogit.mm.map <- function(y, X.re, X.fe, n, shape=1, rate=1,
                       m.0=array(0, dim=c(ncol(X.fe))),
                       P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                       abphi.0=NULL, calc.V=TRUE, trace=FALSE, maxit=100)
{ 
  X = cbind(X.re, X.fe);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(abphi.0)) { abphi.0 = matrix(0, nrow=P+1); abphi.0[P+1] = 1; }

  optim.out <- optim(abphi.0, blogit.llh.mm, gr=NULL, 
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate,
                     m.0=m.0, P.0=P.0,
                     hessian=FALSE,  method="CG",
                     control=list(fnscale=-1, trace=trace, maxit=maxiter));
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

  if (optim.out$convergence != 0) cat("Did not converge:", optim.out$conv, "\n");

  out = list("optim.out"=optim.out, "V"=V, "convergence"=optim.out$convergence, "m"=abphi.pm)
  out$H = hess
  ## out$H.pen = hess.pen

  out
}

blogit.mm.map.3 <- function(y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                            abtm.0=NULL, calc.V=TRUE, trace=FALSE, maxit=1000)
{ 
  X = cbind(X.re, X.fe);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(abtm.0)) { abtm.0 = matrix(0, nrow=P+1+as.numeric(kappa!=0)); }

  optim.out <- optim(abtm.0, blogit.llh.mm.3, hessian=FALSE, method="CG",
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate, kappa=kappa,
                     m.0=m.0, P.0=P.0, control=list(fnscale=-1, trace=trace, maxit=maxit));

  abtm.pm = optim.out$par;
  ## hess = optim.out$hessian
  hess.pen = hessian.blogit.mm.3(abtm.pm, y, X.re, X.fe, n, shape, rate, kappa, P.0)
  hess = hess.pen
  
  V = NA;
  if (calc.V) V = solve(-1 * hess);

  if (optim.out$convergence != 0) cat("Did not converge:", optim.out$conv, "\n");

  out = list("optim.out"=optim.out, "V"=V, "convergence"=optim.out$convergence, "m"=abtm.pm)
  out$H = hess
  ## out$H.pen = hess.pen

  out
}

blogit.mm.map.4 <- function(y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                              m.0=array(0, dim=c(ncol(X.fe))),
                              P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                              abphim.0=NULL, calc.V=TRUE, trace=FALSE, maxit=1000)
{ 
  X = cbind(X.re, X.fe);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(abphim.0)) { abphim.0 = matrix(0, nrow=P+1+as.numeric(kappa!=0)); }
  abphim.0[P+1] = 1;

  optim.out <- optim(abphim.0, blogit.llh.mm.4, gr=grad.blogit.mm.4, 
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate,
                     kappa=kappa, m.0=m.0, P.0=P.0,
                     hessian=FALSE,  method="CG",
                     control=list(fnscale=-1, trace=trace, maxit=maxit));

  ## grad.blogit.mm.4(optim.out$par, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0)
  
  abphim.pm = optim.out$par;
  ## hess = optim.out$hessian
  hess.pen <- hessian.blogit.mm.4(abphim.pm, y, X.re, X.fe, n, shape, rate, kappa, P.0)
  hess = hess.pen

  V = NA;
  if (calc.V) V = solve(-1 * hess);

  if (optim.out$convergence != 0) cat("Did not converge:", optim.out$conv, "\n");

  out = list("optim.out"=optim.out, "V"=V, "convergence"=optim.out$convergence, "m"=abphim.pm)
  out$H = hess
  ## out$H.pen = hess.pen

  out
}

blogit.mm.map.5 <- function(y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                              m.0=array(0, dim=c(ncol(X.fe))),
                              P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                              abtm.0=NULL, calc.V=TRUE, trace=FALSE, maxit=1000)
{
  ## This should be the same as blogt.3 when shape.5 = shape.3 + 2
  X = cbind(X.re, X.fe);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(abtm.0)) { abtm.0 = matrix(0, nrow=P+1+as.numeric(kappa!=0)); }

  optim.out <- optim(abtm.0, blogit.llh.mm.5, hessian=FALSE, method="Nelder-Mead",
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate, kappa=kappa,
                     m.0=m.0, P.0=P.0, control=list(fnscale=-1, trace=trace, maxit=maxit));

  abtm.pm = optim.out$par;
  ## hess = optim.out$hessian
  hess.pen = hessian.blogit.mm.5(abtm.pm, y, X.re, X.fe, n, shape, rate, kappa, P.0)
  hess = hess.pen

  V = NA;
  if (calc.V) V = solve(-1 * hess);

  if (optim.out$convergence != 0) cat("Did not converge:", optim.out$conv, "\n");

  out = list("optim.out"=optim.out, "V"=V, "H"=hess, "H.pen"=hess.pen, "convergence"=optim.out$convergence, "m"=abtm.pm)

}

blogit.mm.map.6 <- function(y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                            dbmt.0=NULL, calc.V=TRUE, trace=FALSE, maxit=1000)
{ 
  X = cbind(X.re, X.fe, 1);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(dbmt.0)) { dbmt.0 = matrix(0, nrow=P+1); }

  optim.out <- optim(dbmt.0, blogit.llh.mm.6, hessian=FALSE, method="CG",
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate, kappa=kappa,
                     m.0=m.0, P.0=P.0, control=list(fnscale=-1, trace=trace, maxit=maxit));

  dbtm.pm = optim.out$par;
  ## hess = optim.out$hessian
  hess.pen = hessian.blogit.mm.6(dbtm.pm, y, X.re, X.fe, n, shape, rate, kappa, P.0)
  hess = hess.pen
  
  V = NA;
  if (calc.V) V = solve(-1 * hess);

  if (optim.out$convergence != 0) cat("Did not converge:", optim.out$conv, "\n");

  out = list("optim.out"=optim.out, "V"=V, "convergence"=optim.out$convergence, "m"=dbtm.pm)
  out$H = hess
  ## out$H.pen = hess.pen

  out
}

blogit.mm.map.7 <- function(y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                              m.0=array(0, dim=c(ncol(X.fe))),
                              P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                              dbmp.0=NULL, calc.V=TRUE, trace=FALSE, maxit=1000)
{ 
  X = cbind(X.re, X.fe, 1);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(dbmp.0)) { dbmp.0 = matrix(0, nrow=P+1); }
  dbmp.0[P+1] = 1;

  optim.out <- optim(dbmp.0, blogit.llh.mm.7, gr=NULL,
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate,
                     kappa=kappa, m.0=m.0, P.0=P.0,
                     hessian=FALSE,  method="CG",
                     control=list(fnscale=-1, trace=trace, maxit=maxit));

  ## grad.blogit.mm.4(optim.out$par, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0)
  
  dbmp.pm = optim.out$par;
  ## hess = optim.out$hessian
  hess.pen <- hessian.blogit.mm.7(dbmp.pm, y, X.re, X.fe, n, shape, rate, kappa, P.0)
  hess = hess.pen

  V = NA;
  if (calc.V) V = solve(-1 * hess);

  if (optim.out$convergence != 0) cat("Did not converge:", optim.out$conv, "\n");

  out = list("optim.out"=optim.out, "V"=V, "convergence"=optim.out$convergence, "m"=dbmp.pm)
  out$H = hess
  ## out$H.pen = hess.pen

  out
}

blogit.mm.map.8 <- function(y, X.re, X.fe, n, shape=1, rate=1, kappa=1,
                            m.0=array(0, dim=c(ncol(X.fe))),
                            P.0=array(0, dim=c(ncol(X.fe), ncol(X.fe))),
                            ab.0=NULL, calc.V=FALSE, trace=FALSE, maxit=1000)
{ 
  X = cbind(X.re, X.fe);
  
  N = nrow(X);
  P = ncol(X);

  if (is.null(ab.0)) { ab.0 = matrix(0, nrow=P); }

  optim.out <- optim(ab.0, blogit.llh.mm.8, gr=NULL,
                     y=y, X.re=X.re, X.fe=X.fe, n=n, shape=shape, rate=rate,
                     kappa=kappa, m.0=m.0, P.0=P.0,
                     hessian=TRUE,  method="CG",
                     control=list(fnscale=-1, trace=trace, maxit=maxit));

  ab.pm = optim.out$par;
  hess  = optim.out$hessian

  V = NA;
  if (calc.V) V = solve(-1 * hess);

  if (optim.out$convergence != 0) cat("Did not converge:", optim.out$conv, "\n");

  out = list("optim.out"=optim.out, "V"=V, "convergence"=optim.out$convergence, "m"=ab.pm)
  out$H = hess
  ## out$H.pen = hess.pen

  out
}

newton.4 <- function(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, abpm.0=NULL,
                     maxiter=100, reltol=1e-8, trace=FALSE)
{
  X = cbind(X.re, X.fe);
  N = nrow(X);
  P = ncol(X);
  P.all = P+1+as.numeric(kappa!=0)

  if (is.null(abpm.0)) { abpm.0 = matrix(0, nrow=P.all); abpm.0[P+1] = 1; }
  if (length(abpm.0) != P.all) return(NA)
  abpm = as.matrix(abpm.0)

  llh.old = blogit.llh.mm.4(abphim=abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0);
  
  go = TRUE
  iter = 0
  
  while (go && iter < maxiter) {
    grad <- grad.blogit.mm.4(abphim=abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0);
    hess <- hessian.blogit.mm.4(abphim=abpm, y, X.re, X.fe, n, shape, rate, kappa, P.0);
    ## Cholesky
    U = chol(-1 * hess);
    V = chol2inv(U)
    ## Eigen
    ## E = eigen(-1 * hess);
    ## V = E$vectors %*% diag(1/E$values) %*% t(E$vectors)
    h = V %*% grad
    abpm = abpm + h
    llh = blogit.llh.mm.4(abphim=abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0);
    if (abs(llh - llh.old) < abs(llh.old) * reltol) go = FALSE
    llh.old = llh
    if (trace != 0) {
      if (iter %% trace == 0) cat("iter:", iter, "llh:", llh, "\n");
    }
    iter = iter + 1
  }

  ## Give things names.
  the.names = c(colnames(X.re), colnames(X.fe), "phi");
  if (kappa != 0) the.names[P+2] = "m"
  abpm = array(abpm, dim=length(abpm));
  names(abpm) = the.names
  
  c.error = 0;
  if (iter >= maxiter) {cat("newton.4: reached max iter.\n"); c.error=1;}

  out = list("hess"=hess, "grad"=grad, "iter"=iter, "c.error"=c.error, "llh"=llh, "mode"=abpm)

  out
}

newton.4.gibbs <- function(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, abpm.0=NULL,
                           maxiter=100, reltol=1e-8, trace=FALSE)
{
  N    = length(y)
  X    = cbind(X.re, X.fe);
  P.a  = ncol(X.re)
  P.b  = ncol(X.fe)
  P.ab = P.a + P.b
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  f.idc = P.ab + 1;
  m.idc = f.idc + 1;
  P.all = P.a + P.b + 1 + as.numeric(kappa!=0)

  if (is.null(abpm.0)) { abpm.0 = matrix(0, nrow=P.ab+1+as.numeric(kappa!=0)); abpm.0[P.ab+1] = 1; }
  if (length(abpm.0) != P.all) return(NA)
  
  abpm = as.numeric(abpm.0)

  llh.old = blogit.llh.mm.4(abphim=abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0);
  
  go = TRUE
  iter = 0
  m  = 4
  maxiter = maxiter * m
  idc = a.idc
  scale = 1
  ## abpm[f.idc] = 1.6 ## checking when phi is fixed.
  
  while (go && iter < maxiter) {
    grad <- grad.blogit.mm.4(abphim=abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0);
    hess <- hessian.blogit.mm.4(abphim=abpm, y, X.re, X.fe, n, shape, rate, kappa, P.0);
    if (iter %% m == 0) {
      idc = a.idc
      scale = 1
    } else if (iter %% m == 1) {
      idc = b.idc
      scale = 1
    } else if (iter %% m == 3) {
      idc = f.idc
      scale = 1
    } else if (iter %% m == 2) {
      idc = m.idc
      scale = 1
    }
    grad.sub = grad[idc]
    hess.sub = hess[idc, idc]
    U = chol(-1 * hess.sub);
    V = chol2inv(U)
    h = V %*% grad.sub
    abpm[idc] = abpm[idc] + h * scale

    llh = blogit.llh.mm.4(abphim=abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0);
    if (abs(llh - llh.old) < abs(llh.old) * reltol) go = FALSE
    llh.old = llh
    if (trace != 0) {
      if (iter %% trace == 0) cat("iter:", iter, "llh:", llh, "\n");
    }
    iter = iter + 1
  }

  ## Give things names.
  the.names = c(colnames(X.re), colnames(X.fe), "phi");
  if (kappa != 0) the.names[m.idc] = "m"
  abpm = array(abpm, dim=P.all);
  names(abpm) = the.names

  c.error = 0;
  if (iter >= maxiter) {cat("newton.4.gibbs: reached max iter.\n"); c.error=1;}

  out = list("hess"=hess, "grad"=grad, "iter"=iter, "c.errer"=c.error, "llh"=llh, "mode"=abpm)

  out
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

  out <- list(abphi = array(0, dim=c(samp, P.ab+1)),
              a.prob = rep(0, samp)
              )

  ## Find mode
  out.map = blogit.mm.map(y, X.re, X.fe, n, shape, rate, m.0, P.0, abphi.0=NULL, maxit=100000)
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
      out$abphi[i-burn,]    = out.abphi$abphi
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
  out$optim      = out.map

  out$ab  = out$abphi[,1:P.ab]
  out$phi = out$abphi[,P.ab+1]

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

  out <- list(abphi = array(0, dim=c(samp, P.ab+1)),
              a.prob = rep(0, samp)
              )

  ## Find mode
  out.map = blogit.mm.map(y, X.re, X.fe, n, shape, rate, m.0, P.0, abphi.0=NULL, maxit=100000)
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
      out$abphi[i-burn,]    = out.abphi$abphi
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

  out$ab  = out$abphi[,1:P.ab]
  out$phi = out$abphi[,P.ab+1]

  out
}

################################################################################
                              ## IND. METROP 3 ##
################################################################################

log.tget.to.ppsl.blogit.3 <- function(abtm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, m, UorP, df=Inf, is.prec=TRUE)
{
  ## Calculate log f - log q
  if (df==Inf)
    log.fdivq = llh.norm(abtm, m, UorP, is.prec)
  else
    log.fdivq = llh.t(abtm, m, UorP, df, is.prec)

  log.fdivq = blogit.llh.mm.3(abtm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0) - log.fdivq
  
  log.fdivq
}

draw.abtm.ind.MH.3 <- function(abtm, y, X.re, X.fe, n, shape, rate, kappa, m0, P0, map, U.map, log.fdivq, df=Inf)
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
  
  P = length(abtm)
  zero = rep(0, P);

  ## log.fdivq = log.tget.to.ppsl(beta, y, X, n, m.0, P.0, map, U.map, df=df, is.prec=FALSE)

  ## Propose
  if (df==Inf) ep.ppsl = rnorm(P) else ep.ppsl = rt(P, df=df)

  ppsl = map + backsolve(U.map, ep.ppsl);

  log.fdivq.ppsl <- log.tget.to.ppsl.blogit.3(ppsl, y, X.re, X.fe, n, shape, rate, kappa,
                                              m0, P0, map, U.map, df=df, is.prec=FALSE)

  ## acceptance prob
  a.prob = min( exp(log.fdivq.ppsl - log.fdivq), 1);
  accept = runif(1) < a.prob

  if (accept) {
    abtm = ppsl;
    log.fdivq = log.fdivq.ppsl
  }

  out = list("abtm"=abtm, "a.prob"=a.prob, "accept"=accept, "log.fdivq"=log.fdivq)

  out
}

ind.metropolis.blogit.3 <- function(y, X.re, X.fe, n=rep(1, length(y)), shape=1, rate=1, kappa=1,
                                    m.0=rep(0, ncol(X.fe)), P.0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe)),
                                    samp=1000, burn=100, verbose=1000, df=Inf)
{
  ## An independent MH routine that can be adapted to varying P0.
  ## Assume proper posterior
  
  ## y: response, number of responses in each category
  ## X: design
  ## n: number of trials per draw.
  ## kappa: precision inflation.
  
  y = as.matrix(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  
  N = nrow(X.re);
  P.a  = ncol(X.re); 
  P.b  = ncol(X.fe); 
  P.ab = P.a + P.b
  P.all = P.ab + 1 + as.numeric(kappa!=0);

  a.idc = 1:P.a;
  b.idc = 1:P.b + P.a;
  t.idc = P.ab + 1
  m.idc = P.ab + 2

  ## n = rep(1, N)
  
  if (is.null(m.0)) m.0 = rep(0, P.b);
  if (is.null(P.0)) {
    p.0 = 0e-4
    P.0 = diag(p.0, P.b);
  }

  out <- list(abtm = matrix(0, nrow=samp, ncol=P.all),
              a.prob = rep(0, samp)
              )

  ## Find mode
  out.map = blogit.mm.map.3(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, abtm.0=NULL, maxit=100000)
  m = out.map$m
  U = chol(-1 * out.map$H)
  cat("Finished with MLE.  Convergence =", out.map$convergence, "\n")

  ## start at ppsl mean.
  ep   = rep(0, P.all)
  abtm = out.map$m
  log.fdivq = log.tget.to.ppsl.blogit.3(abtm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, m, U, df=df, is.prec=FALSE)
  
  ## Generate proposal
  out.abtm = draw.abtm.ind.MH.3(abtm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, map=m, U.map=U, log.fdivq=log.fdivq, df=df)
  abtm = out.abtm$abtm
  log.fdivq = out.abtm$log.fdivq

  ## Timing
  start.time = proc.time()
  naccept = 0

  ## Do Metropolis ##
  for (i in 1:(samp+burn)) {

    if (i==burn+1) { start.ess = proc.time(); }

    ## Generate proposal
    out.abtm = draw.abtm.ind.MH.3(abtm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, map=m, U.map=U, log.fdivq=log.fdivq, df=df)
    abtm = out.abtm$abtm
    log.fdivq = out.abtm$log.fdivq
    
    if (i > burn) {
      out$abtm[i-burn, ] = abtm
      out$a.prob[i-burn] = out.abtm$a.prob
      naccept = naccept + out.abtm$accept
    }

    if (i %% verbose == 0) {
      if (i > burn) cat("Ind MM3: Ave a.prob:", mean(out$a.prob[1:(i-burn)]), ", ");
      cat("Iteration:", i, "\n");
    }
    
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$naccept    = naccept
  out$acceptr    = naccept / samp;
  out$optim      = out.map

  out$abm   = out$abtm[,c(a.idc, b.idc, ifelse(kappa!=0, m.idc, NULL))]
  out$phi   = exp(out$abtm[,t.idc]);

  out
}

################################################################################
                              ## IND. METROP 4 ##
################################################################################

log.tget.to.ppsl.blogit.4 <- function(abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, m, UorP, df=Inf, is.prec=TRUE)
{
  ## Calculate log f - log q
  if (df==Inf)
    log.fdivq = llh.norm(abpm, m, UorP, is.prec)
  else
    log.fdivq = llh.t(abpm, m, UorP, df, is.prec)

  log.fdivq = blogit.llh.mm.4(abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0) - log.fdivq
  
  log.fdivq
}

draw.abpm.ind.MH.4 <- function(abpm, y, X.re, X.fe, n, shape, rate, kappa, m0, P0, map, U.map, log.fdivq, df=Inf)
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
  
  P = length(abpm)
  zero = rep(0, P);

  ## log.fdivq = log.tget.to.ppsl(beta, y, X, n, m.0, P.0, map, U.map, df=df, is.prec=FALSE)

  ## Propose
  if (df==Inf) ep.ppsl = rnorm(P) else ep.ppsl = rt(P, df=df)

  ppsl = map + backsolve(U.map, ep.ppsl);

  log.fdivq.ppsl <- log.tget.to.ppsl.blogit.4(ppsl, y, X.re, X.fe, n, shape, rate, kappa,
                                              m0, P0, map, U.map, df=df, is.prec=FALSE)

  ## acceptance prob
  a.prob = min( exp(log.fdivq.ppsl - log.fdivq), 1);
  accept = runif(1) < a.prob

  if (accept) {
    abpm = ppsl;
    log.fdivq = log.fdivq.ppsl
  }

  out = list("abpm"=abpm, "a.prob"=a.prob, "accept"=accept, "log.fdivq"=log.fdivq)

  out
}

ind.metropolis.blogit.4 <- function(y, X.re, X.fe, n=rep(1, length(y)), shape=1, rate=1, kappa=1,
                                    m.0=rep(0, ncol(X.fe)), P.0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe)),
                                    samp=1000, burn=100, verbose=1000, df=Inf, thin=1, center=NULL)
{
  ## An independent MH routine that can be adapted to varying P0.
  ## Assume proper posterior
  
  ## y: response, number of responses in each category
  ## X: design
  ## n: number of trials per draw.
  ## kappa: precision inflation.
  thin = round(thin)
  if (thin<1) {cat("thin must be > 0.\n"); return(NA); }
  
  y = as.matrix(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  
  N = nrow(X.re);
  P.a  = ncol(X.re); 
  P.b  = ncol(X.fe); 
  P.ab = P.a + P.b
  P.all = P.ab + 1 + as.numeric(kappa!=0);

  a.idc = 1:P.a;
  b.idc = 1:P.b + P.a;
  t.idc = P.ab + 1
  m.idc = P.ab + 2

  ## n = rep(1, N)
  
  if (is.null(m.0)) m.0 = rep(0, P.b);
  if (is.null(P.0)) {
    p.0 = 0e-4
    P.0 = diag(p.0, P.b);
  }

  out <- list(abpm = matrix(0, nrow=samp, ncol=P.all),
              a.prob = rep(0, samp)
              )
  
  ## Set the proposal.
  if (!is.null(center)) {
    grad <- grad.blogit.mm.4(abphim=center, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0);
    hess <- hessian.blogit.mm.4(abphim=center, y, X.re, X.fe, n, shape, rate, kappa, P.0);
    ## Cholesky
    U = chol(-1 * hess);
    V = chol2inv(U)
    m = V %*% grad
    out.map = list("c"=0)
  } else {
    ## Find mode
    ## out.map = blogit.mm.map.4(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, abphim.0=NULL, maxit=100000)
    ## m = out.map$m
    ## U = chol(-1 * out.map$H)
    out.map <- newton.4.gibbs(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, abpm.0=NULL,
                              maxiter=100000, reltol=1e-14, trace=FALSE)
    m = out.map$m
    U = chol(-1 * out.map$hess);
  }

  cat("Finished with MLE.  Convergence =", out.map$c, "\n")

  ## start at ppsl mean.
  ep   = rep(0, P.all)
  abpm = m
  log.fdivq = log.tget.to.ppsl.blogit.4(abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, m, U, df=df, is.prec=FALSE)
  
  ## Generate proposal
  out.abpm = draw.abpm.ind.MH.4(abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, map=m, U.map=U, log.fdivq=log.fdivq, df=df)
  abpm = out.abpm$abpm
  log.fdivq = out.abpm$log.fdivq

  ## Total samp/burn
  tsamp = thin * samp
  tburn = thin * burn
  tverb = thin * verbose
  
  ## Timing
  start.time = proc.time()
  naccept = 0

  ## Do Metropolis ##
  for (i in 1:(tsamp+tburn)) {

    if (i==tburn+1) { start.ess = proc.time(); }

    ## Generate proposal
    out.abpm = draw.abpm.ind.MH.4(abpm, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, map=m, U.map=U, log.fdivq=log.fdivq, df=df)
    abpm = out.abpm$abpm
    log.fdivq = out.abpm$log.fdivq
    
    if (i > tburn && i %% thin == 0) {
      out$abpm[i/thin-burn, ] = abpm
      out$a.prob[i/thin-burn] = out.abpm$a.prob
    }
    naccept = naccept + out.abpm$accept * as.numeric(i > tburn)

    if (i %% tverb == 0) {
      if (i > tburn) cat("Ind MM4: a.rate: ", naccept, "/", (i-tburn),
                         " = ", signif(naccept/(i-tburn), 4), ", ", sep="");
      cat("Iteration:", i / thin, "\n");
    }
    
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$naccept    = naccept / thin
  out$acceptr    = naccept / tsamp;
  out$optim      = out.map

  out$phi   = out$abpm[,t.idc]
  out$abm   = out$abpm[,c(a.idc, b.idc, ifelse(kappa!=0, m.idc, NULL))]

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
        output$w[j-burn,]  = w
        output$ab[j-burn,] = ab
        output$phi[j-burn] = phi
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
} ## logit.PG.mm

################################################################################
                              ## Logit PG MM 2 ##
################################################################################

## Bayesian logistic regression
##------------------------------------------------------------------------------
logit.PG.mm.2 <- function(y, X.re, X.fe, n=rep(1, length(y)), shape=1, rate=1, kappa=1,
                        m0=rep(0, ncol(X.fe)), P0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe)),
                        samp=1000, burn=500, verbose=500, seed=list("abm"=NULL, "phi"=NULL))
{
  ## X: n by p matrix
  ## y: n by 1 vector, total # response
  ## n: n by 1 vector, # of obs at distinct x

  y = as.numeric(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  X    = cbind(X.re, X.fe)

  ## Adjust if we are using m.
  inc.m = as.numeric(kappa != 0);
  if (inc.m) X = cbind(X, 1)
  ## When kappa == Inf, we need to take out phi^1/2.  kappa=Inf will take care of m^2.
  not.flat = as.numeric(kappa != Inf);
  
  P.a   = ncol(X.re)
  P.b   = ncol(X.fe)
  P.ab  = P.a + P.b
  P.all = P.ab + inc.m
  N = nrow(X)
  
  a.idc = 1:P.a
  b.idc = 1:P.b + P.a
  m.idc = P.ab + 1

  ## Precompute
  Z = colSums(X * (y-n/2));
  Z[b.idc] = Z[b.idc] + P0 %*% m0;
  s.post = shape + P.a + inc.m * not.flat

  ## PsiToBeta = solve(t(X) %*% X) %*% t(X);

  w     = rep(0,N)
  ## w = w.known;
  dbm   = rep(0.0, P.all)
  m     = 0
  phi   = shape / rate;

  if (!is.null(seed$abm)) { dbm = seed$abm; dbm[a.idc] = seed$abm[a.idc] - seed$abm[m.idc]; }
  if (!is.null(seed$phi)) { phi = seed$phi; }

  out <- list(w = matrix(nrow=samp, ncol=N),
              abm = matrix(nrow=samp, ncol=P.ab+1),
              phi = rep(0, samp)
              )

  ## Timing
  start.time = proc.time()
  
  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();
    
    ## draw w
    psi = drop(X %*% dbm)
    
    ## Devroye is faster anyway.
    w = rpg.devroye(N, n, psi);
    
    ## Draw (alpha, beta, m) - Joint Sample.
    ## Draw (delta = alpha - m, beta, m), then adjust.
    PP = t(X) %*% (X * w);
    PP[b.idc, b.idc] = PP[b.idc,b.idc] + P0
    PP[a.idc, a.idc] = PP[a.idc, a.idc] + diag(phi, P.a);
    if (kappa != 0) PP[m.idc, m.idc] = PP[m.idc,m.idc] + phi / kappa^2;
    ## U = chol(PP);
    ## pm = backsolve(U, Z, transpose=TRUE);
    ## pm = backsolve(U, m);
    ## abm = pm + backsolve(U, rnorm(P.ab))
    S = chol2inv(chol(PP));
    pm = S %*% as.vector(Z);
    dbm = pm + t(chol(S)) %*% rnorm(P.all);

    delta = dbm[a.idc]
    if (inc.m) m = dbm[m.idc];
    abm = dbm; abm[a.idc] = delta + m;

    ## Draw phi
    r.post = t(delta) %*% (delta) + rate + ifelse(inc.m, m^2 / kappa^2, 0)
    phi = rgamma(1, shape=0.5*s.post, rate=0.5*r.post);
    
    # Record if we are past burn-in.
    if (j>burn) {
        out$w[j-burn,]   = w
        out$abm[j-burn,] = abm
        out$phi[j-burn]  = phi
    }

    if (j %% verbose == 0) { print(paste("LogitPG MM 2: Iteration", j)); }
  }

  if (inc.m) out$fe = out$abm[,c(m.idc, b.idc)]
  else out$fe = cbind(0, out$abm[,b.idc])
  out$re = out$abm[,a.idc] - out$fe[,1]
  out$delta = out$re
  out$alpha = out$abm[,a.idc]

  colnames(out$fe) = c("Intercept", colnames(X.fe))
  
  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess

  ## ## Add new data to out.
  ## out$"y" = y;
  ## out$"X" = X;
  ## out$"n" = n;

  out
} ## logit.PG.mm.2

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
  
  out.ind.3 = ind.metropolis.blogit.3(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, samp=samp, burn=burn, verbose=verbose, df=df)
  out.ind.4 = ind.metropolis.blogit.4(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, samp=samp, burn=burn, verbose=verbose, df=df)
  out.pg.mm.2 = logit.PG.mm.2(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, samp=samp, burn=burn, verbose=verbose);

  out.mlogit = mlogit.MH.R(y, X, n, m.0=m.0.mlogit, P.0=P.0.mlogit, beta.0=NULL,
    samp=samp, burn=burn, method="Ind", tune=1.0, df=df, verbose=1000)
  
  out.pg = logit.R(y, X, m0=m.0.PG.mm, P0=P.0.PG.mm, samp=samp, burn=burn, verbose=verbose);
  
  rbind(sum.stat(out.ind$ab, out.ind$ess.time[3]),
        sum.stat(out.ind$phi, out.ind$ess.time[3]))

  rbind(sum.stat(out.ind.2$ab, out.ind.2$ess.time[3]),
        sum.stat(out.ind.2$phi, out.ind.2$ess.time[3]))

  rbind(rbind(sum.stat(out.pg.mm$ab, out.pg.mm$ess.time[3])),
        rbind(sum.stat(out.pg.mm$phi, out.pg.mm$ess.time[3])))
  
  rbind(sum.stat(out.ind.3$abm, out.ind.3$ess.time[3]),
        sum.stat(out.ind.3$phi, out.ind.3$ess.time[3]))
  
  rbind(sum.stat(out.ind.4$abm, out.ind.4$ess.time[3]),
        sum.stat(out.ind.4$phi, out.ind.4$ess.time[3]))

  rbind(rbind(sum.stat(out.pg.mm.2$abm, out.pg.mm$ess.time[3])),
        rbind(sum.stat(out.pg.mm.2$phi, out.pg.mm$ess.time[3])))

  rbind(sum.stat(out.pg$beta, out.pg$ess.time[3]))

  rbind(sum.stat(out.mlogit$beta, out.mlogit$ess.time[3]))

  sstat.ind.1 = sum.stat(out.ind$ab, out.ind$ess.time[3]);
  sstat.ind.2 = sum.stat(out.ind.2$ab, out.ind$ess.time[3]);
  sstat.pg.mm = sum.stat(out.pg.mm$ab, out.pg.mm$ess.time[3]);
  
}

################################################################################
                           ## MIXED MODEL EXAMPLE ##
################################################################################

if (FALSE)
{

  library(mlmRev)
  library(lme4)
  data(Contraception)
  summary(Contraception)
  fm1 <- glmer(use ~ urban+age+livch+(1|district), Contraception, binomial)
  summary(fm1)
  ranef(fm1)

  y    = as.numeric(Contraception$use=="Y")
  X.fe = model.matrix(use ~ urban + age + livch, data=Contraception)
  X.fe = X.fe[,-1]
  ## District 54 is missing.
  X.re = model.matrix(use ~ factor(district) + 0, data=Contraception)
  ## Add it back?
  ## X.re = cbind(X.re[,1:35], 0, X.re[,55:61])
  X = cbind(X.re, X.fe)
  
  N = nrow(X.fe)
  P = ncol(X.fe)
  n = rep(1, N)
  P.PG = ncol(X)
  
  shape = 2
  rate =  2
  kappa = Inf
  inc.m = as.numeric(kappa!=0)
  m.0 = rep(0, P)
  P.0 = diag(0.01, P)
  
  m0.PG = rep(0, P.PG)
  P0.PG = diag(0.01, P.PG)
  
  samp = 2000
  burn = 2000
  verbose = 1000

  ## source("Logit-MixedModel.R")
  out.0 = logit.R(y, X, m0=m0.PG, P0=P0.PG, samp=samp, burn=burn, verbose=verbose);
  out.1 = logit.PG.mm(y, X.re, X.fe, n, shape, rate, m.0, P.0, samp=samp, burn=burn, verbose=verbose);
  out.2 = logit.PG.mm.2(y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0, samp=samp, burn=burn, verbose=verbose);
  colnames(out.2$re) = paste("D", c(1:53, 55:61), sep="")
  
  colMeans(out.2$re)
  ranef(fm1)$district[,1]
  colMeans(out.2$fe)
  summary(fm1)

  ## png(filename="ranef.png", width=800, height=200)
  par(mfrow=c(1,1))
  par(mar=c(4,4,3,3))
  boxplot(out.0$beta[,1:60], xlab="District", ylab="Random Effect", main="Random Intercept by District", names=colnames(out.2$re))
  dev.off()
  
  ## png(filename="mm-ranef.png", width=800, height=200)
  par(mfrow=c(1,1))
  par(mar=c(4,4,3,3))
  boxplot(out.2$alpha, xlab="District", ylab="Random Int.", main="Random Intercept by District", names=colnames(out.2$re))
  dev.off()

  id = matrix(c(1:53, 55:61), nrow=samp, ncol=60, byrow=TRUE)
  plot(id, out.2$re, col="#11111111", pch=19)

  par(mfrow=c(2,3))
  for(i in 1:ncol(out.2$fe)) {
    the.name = colnames(out.2$fe)[i]
    hist(out.2$fe[,i], main=the.name, ylab="", xlab=the.name, prob=TRUE)
  }

  
  
}
