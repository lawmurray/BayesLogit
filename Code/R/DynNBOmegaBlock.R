
source("DynLogitBetaBlockMH.R")
source("NB-Shape.R")

dyn.nb.om <- function(y, X.dyn, m0, C0,
                         samp=1000, burn=100, verbose=100, starts=c(1),
                         mu.m0=NULL, mu.P0=NULL,
                         phi.m0=NULL, phi.P0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         X.stc=NULL, m0.stc=NULL, C0.stc=NULL,
                         mu.true = NULL, phi.true=NULL, W.true=NULL,
                         alpha.true=NULL, beta.true=NULL, d.true=NULL,
                         just.max=FALSE)
{
  ## Dim and restructure.
  X.dyn = as.matrix(X.dyn);
    
  C0  = as.matrix(C0);
  T   = length(y)
  N.b = ncol(X.dyn)
  M   = samp
  n.starts = length(starts)

  if (is.null(X.stc)) {
    X.stc      = matrix(1, T, 1);
    alpha.true = 0;
  }
  N.a = ncol(X.stc);

  tX.dyn = t(X.dyn);
  
  ## Default prior parameters -- almost a random walk for beta ##
  if (is.null(m0)     || is.null(C0))     { m0     = rep(0.0, N.b) ; C0     = diag(1.0, N.b); }
  if (is.null(mu.m0)  || is.null(mu.P0))  { mu.m0  = rep(0.0 ,N.b) ; mu.P0  = rep(100 , N.b); }
  if (is.null(phi.m0) || is.null(phi.P0)) { phi.m0 = rep(0.99,N.b) ; phi.P0 = rep(100 , N.b); }
  if (is.null(W.a0)   || is.null(W.b0))   { W.a0   = rep(1.0, N.b) ; W.b0   = rep(1.0,  N.b);  }
  if (is.null(m0.stc) || is.null(C0.stc)) { m0.stc = rep(0  , N.a) ; C0.stc = diag(100, N.a); }
  
  ## Output
  out <- list("beta"  = array(0, dim=c(M, N.b, T)),
              "omega" = array(0, dim=c(M, N.b, T)),
              "alpha" = array(0, dim=c(M, N.a)),
              "psi"   = array(0, dim=c(M, T)),
              "mu"    = array(0, dim=c(M, N.b)),
              "phi"   = array(0, dim=c(M, N.b)),
              "W"     = array(0, dim=c(M,N.b)),
              ## "ppsl.b"=array(0, dim=c(M, N.b, T+1)),
              "ac.rate" = array(0, dim=c(M, 2)))

  ## Preprocess ## 
  ymax = max(y);
  F = cumsum(hist(y, breaks=0:(ymax+1)-0.5, plot=FALSE)$counts)
  G = T - F;
  
  ## Initialize
  beta  = array(0, dim=c(N.b, T));
  alpha = array(0, dim=c(N.a));
  omega = beta;
  d     = 1

  mu   = mu.m0
  phi  = phi.m0
  W    = W.b0 / W.a0;

  naccept = c(0,0)

  P0.stc = solve(C0.stc);
  b0.stc = P0.stc %*% m0.stc;

  psi.dyn = colSums(tX.dyn * beta);
  psi.stc = X.stc %*% alpha;

  ## Check if known.
  know.phi <- know.mu <- know.W <- know.d <- know.alpha <- know.beta <- FALSE
  
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (any(phi==1)) {  mu.true   = rep(0, N.b); } }
  if (!is.null(mu.true))   { mu   = mu.true ;  know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true  ;  know.W    = TRUE; }
  if (!is.null(d.true))    { d    = d.true  ;  know.d    = TRUE; }
  
  if (!is.null(alpha.true)) {
    alpha = alpha.true;
    psi.stc  = X.stc %*% alpha;
    know.alpha = TRUE;
  }

  if (!is.null(beta.true)) {
    beta = beta.true
    psi.dyn = colSums(tX.dyn * beta);
    know.beta = TRUE
  }

  ## Check
  ## cat("psi.stc", psi.stc, "\n");
    
  ## if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  llh = nb.mu.llh.3(y, psi.dyn, psi.stc, d=d)
  
  start.time = proc.time()
  
  ## MCMC
  for(i in 1:(samp+burn)) {
    if (i==burn+1) {
      start.ess = proc.time();
      naccept = c(0,0);
    }

    ## Draw beta
    if (!know.beta) {
      prior.prec = diag(1/W, N.b * T);
      prior.prec[1:N.b,1:N.b] = solve(diag(W, N.b) + diag(phi,N.b) %*% C0 %*% diag(phi,N.b));
      draw = draw.omega.mh(omega, beta, llh, y, tX.dyn, ntrials=d, prior.prec, phi, starts, just.max=just.max, offset=psi.stc, type=1);
      omega   = draw$omega;
      beta    = draw$beta;
      llh     = draw$llh
      psi.dyn = colSums(tX.dyn * beta)
      naccept[1] = naccept[1] + draw$naccept / n.starts;
    }
    
    ## Draw alpha
    if (!know.alpha) {
      draw = draw.stc.beta.mh(beta=alpha, llh, y, X.stc, ntrials=d, b0=b0.stc, P0=P0.stc, just.max=just.max, offset=psi.dyn, type=1)
      alpha    = draw$beta;
      llh      = draw$llh;
      psi.stc  = X.stc %*% alpha;
      naccept[2] = naccept[2] + draw$naccept
    }

    ## Draw d
    if (!know.d) {
      d = draw.df(y, d, exp(llh$psi), G, ymax);
      ## Need to do this since d changes.
      llh = nb.mu.llh.3(y, psi.dyn, psi.stc, d=d)
    }              
    
    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    ## W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    if (!know.mu)  mu  = draw.mu.ar1.ind (beta, phi, W, mu.m0, mu.P0)
    if (!know.phi) phi = draw.phi.ar1.ind(beta, mu, W, phi.m0, phi.P0, phi)
    if (!know.W)   W   = draw.W.ar1.ind  (beta, mu, phi, W.a0, W.b0)
    
    if (i > burn) {
      ii = i - burn
      ## out$ppsl.b[i-burn,,] = ppsl.b
      out$beta[ii,,]   = beta
      out$omega[ii,,]  = omega
      out$alpha[ii,]   = alpha
      out$psi[ii,]     = llh$psi
      out$mu[ii,]      = mu
      out$phi[ii,]     = phi
      out$W[ii,]       = W;
      out$ac.rate[ii,]  = naccept / ii
    }

    if (i %% verbose == 0) {
      cat("DynNBOmega: ");
      if (i > burn) { cat("Accept rate:", naccept / (i-burn), ", ") }
      cat("Iteration:", i, "\n");
    }
    
  }
              
  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$last       = draw
  out$a.rate     = naccept / (i - burn);
  out$error      = 0
  
  out
  
}

################################################################################

if (FALSE) {

  N = 400
  P = 1

  P.a = 2
  X.stc = cbind(1, rnorm(N, mean=0, sd=0.2));
  alpha.t = c(1, -1);
  ## alpha.t = c(0, 0);
  
  phi = 0.9
  sig = 0.3
  W   = sig^2;
  n = 4
  m0.dyn = rep(0, P)
  C0.dyn = diag(1.0,P)
  m0.stc = rep(0, P.a)
  C0.stc = diag(100, P.a)

  beta.t  = matrix(0, nrow=P, ncol=N);
  omega.t = matrix(rnorm(N*P, 0, sig), nrow=P, ncol=N);
  
  beta.1  = t(chol(C0)) %*% rnorm(P) + rnorm(P, 0, sig);
  beta.t[,1]  = beta.1;
  omega.t[,1] = beta.1;
  
  for (i in 2:N) {
    beta.t[,i] = phi * beta.t[,i-1] + omega.t[,i];
  }

  ## psi.t is the LOG MEAN here.
  X = matrix(1, nrow=N, ncol=P)
  psi.dyn.t = colSums(t(X) * beta.t)
  psi.stc.t = X.stc %*% alpha.t
  psi.t = psi.dyn.t + psi.stc.t;
  mu.t  = exp(psi.t);
  y = rnbinom(N, n, mu=mu.t);

  plot(y)
  lines(mu.t, col=1, lty=2);
  
}

if (FALSE) {

  ## starts = c(1,5)
  starts = c(1,6);
  samp = 10
  burn = 0
  verbose = 1000
  just.max = FALSE

  phi.in = NULL
  W.in   = NULL

  Ny = length(y)
  starts = seq(1, Ny, 2)
  
  ## set.seed(1000);
  
  source("DynNBOmegaBlock.R")
  start.time = proc.time()
  output4 <- dyn.nb.om(y=y, X.dyn=X, m0=m0.dyn, C0=C0.dyn,
                       samp=samp, burn=burn, verbose=verbose, starts=starts,
                       mu.m0=NULL, mu.P0=NULL,
                       phi.m0=NULL, phi.P0=NULL,
                       W.a0=NULL, W.b0=NULL,
                       X.stc=X.stc, m0.stc=m0.stc, C0.stc=C0.stc,
                       mu.true = 0, phi.true=phi.in, W.true=W.in,
                       alpha.true=NULL, beta.true = NULL, d.true=n, 
                       just.max=just.max)
  diff.time = proc.time() - start.time
  print(diff.time)
  
  ## plot(beta.t[1,], type="l")
  ## lines(output4$beta[samp,1,], col=2)

  alpha.m4 = apply(output4$alpha, 2, mean)
  beta.m4  = apply(output4$beta, c(2,3), mean)
  ## psi.m4   = colSums(t(X) * beta.m4) + X.stc %*% alpha.m4
  psi.m4   = output4$psi
  mu.m4    = apply(exp(output4$psi), 2, mean)

  plot(y)
  lines(mu.t, col=1)
  lines(mu.m4, col=4, lty=4)

  ##------------------------------------------------------------------------------
  
  pg.C0 = diag(1, P + P.a); pg.C0[1:P.a,1:P.a] = C0.stc; pg.C0[(P.a+1):(P.a+P),(P.a+1):(P.a+P)] = C0.dyn
  pg.m0 = c(m0.stc, m0.dyn)

  ## pg.C0 = C0
  ## pg.m0 = m0

  source("DynNBPG.R")
  start.time = proc.time()
  out.pg <- dyn.NB.PG(y=y, X.dyn=X, X.stc=X.stc,
                      samp=samp, burn=burn, verbose=verbose,
                      m.0=pg.m0, C.0=pg.C0,
                      mu.m0=NULL, mu.P0=NULL,
                      phi.m0=NULL, phi.P0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      d.true=n, w.true=NULL,
                      beta.true=NULL, iota.true=NULL,
                      mu.true=0, phi.true=phi, W.true=W)
  diff.time = proc.time() - start.time
  print(diff.time)
    
  beta.pg  = apply(out.pg$beta[,,-1,drop=FALSE], c(2,3), mean)
  alpha.pg = apply(out.pg$alpha, 2, mean)
  ## psi.m5 = colSums(t(X) * beta.pg) + X.stc %*% alpha.pg
  mu.m5    = apply(exp(out.pg$lmean), 2, mean)

  lines(mu.m5, col=5, lty=5)
  
}
