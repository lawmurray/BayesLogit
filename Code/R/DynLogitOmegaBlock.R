

dyn.logit.om <- function(y, X.dyn, n, m0, C0,
                         samp=1000, burn=100, verbose=100, starts=c(1),
                         mu.m0=NULL, mu.P0=NULL,
                         phi.m0=NULL, phi.P0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         X.stc=NULL, m0.stc=NULL, C0.stc=NULL,
                         mu.true = NULL, phi.true=NULL, W.true=NULL,
                         alpha.true=NULL, beta.true=NULL,
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

  ## Initialize
  beta  = array(0, dim=c(N.b, T));
  alpha = array(0, dim=c(N.a));
  omega = beta;

  mu   = mu.m0
  phi  = phi.m0
  W    = W.b0 / W.a0;

  naccept = c(0,0)

  P0.stc = solve(C0.stc);
  b0.stc = P0.stc %*% m0.stc;

  psi.dyn = colSums(tX.dyn * beta);
  psi.stc = X.stc %*% alpha;

  llh = logit.llh.3(y, psi.dyn, psi.stc, n=n)

  ## Check if known.
  know.phi <- know.mu <- know.W <- know.alpha <- know.beta <- FALSE
  
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (any(phi==1)) {  mu.true   = rep(0, N.b); } }
  if (!is.null(mu.true))   { mu   = mu.true ;  know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true  ;  know.W    = TRUE; }
  
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

  start.time = proc.time()
  
  ## MCMC
  for(i in 1:(samp+burn)) {
    if (i==burn+1) {
      start.ess = proc.time();
      naccept = c(0,0);
    }
    
    ## Draw beta
    if (!know.beta) {
      prior.prec = diag(1/W, N.b * N);
      prior.prec[1:N.b,1:N.b] = solve(diag(1/W, N.b) + diag(phi,P) %*% C0 %*% diag(phi,P));
      draw = draw.omega.mh(omega, beta, llh, y, tX.dyn, ntrials=n, prior.prec, phi, starts, just.max=just.max, offset=psi.stc);
      omega   = draw$omega;
      beta    = draw$beta;
      llh     = draw$llh
      psi.dyn = colSums(tX.dyn * beta)
      naccept[1] = naccept[1] + draw$naccept / n.starts;
    }
    
    ## Draw alpha
    if (!know.alpha) {
      draw = draw.stc.beta.mh(beta=alpha, llh, y, X.stc, ntrials=n, b0=b0.stc, P0=P0.stc, just.max=just.max, offset=psi.dyn)
      alpha    = draw$beta;
      llh      = draw$llh;
      psi.stc  = X.stc %*% alpha;
      naccept[2] = naccept[2] + draw$naccept
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
      cat("DynLogitOmega: ");
      if (i > burn) { cat("Accept rate:", naccept / (i-burn), ", ") }
      cat("Iteration:", i, "\n");
    }
    
  }
              
  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$last       = draw
  out$a.rate     = naccept / (i - burn);
  
  out
  
}

################################################################################

