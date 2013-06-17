
## llh and llh.2 defined in 
source("DynLogitBetaBlockMH.R")

get.laplace.2 <- function(omega, y, XiL, df, D, block.start, block.end, bsize)
{
  ## XiL is diag(x_t') L
  ## D is diagonal precision.

  nblocks = block.end - block.start + 1;
  inb.idc = block.start:block.end
  ## inm.idc = (bsize * (block.start-1) + 1):(bsize * block.end);
  inm.idc = bidc.to.midc(block.start, bsize, block.end);

  ## cat("inb.idc:", inb.idc, "\n");
  ## cat("inm.idc:", inm.idc, "\n");
  
  Gtil = t(df$l1)  %*% XiL;
  Ftil = t(df$phi) %*% XiL;
  btil = as.numeric(Gtil - Ftil);
  Htil = t(XiL) %*% (XiL * df$l2);
  Otil = D - Htil;

  Om.11 = Otil[inm.idc, inm.idc, drop=FALSE];

  omw.12 = 0.0
  if (length(inm.idc) != nrow(D)) {
    Om.12   = Otil[inm.idc,-inm.idc,drop=FALSE]
    omega.2 = as.numeric( omega[,-inm.idc] );
    omw.12  = Om.12 %*% omega.2;
    ## cat("Om.12:\n"); print(Om.12);
    ## cat("omega.2:", omega.2, "\n");
    ## cat("omw.12:", omw.12, "\n");
  }

  b1 = btil[inm.idc];
  ## cat("b1:", b1, "\n");

  ## cat("df$l1:", df$l1, "\n");
  ## cat("Gtil:", Gtil, "\n");
  ## cat("Ftil:", Ftil, "\n");
  ## cat("omw.12:", omw.12, "\n");
  ## cat("b1:", b1, "\n");

  V.post = solve(Om.11)
  m.post = - V.post %*% (omw.12 - b1);

  ## cat("m:", m.post, "\n"); 
  ## ##print(V.post)
  ## print(Om.11);
  ## print(chol(Om.11));
  
  out = list("P"=Om.11, "V"=V.post, "m"=m.post);
  out 
}

draw.omega.block.2 <- function(y, X, omega, df, XiL, L, DP, block.start, block.end, n=1, just.max=FALSE)
{
  ## Assume X is N x P like a design matrix.
  ## DP is prior precision for omega.
  
  inb.idc = block.start:block.end;
  nblocks = block.end - block.start + 1;
  bsize   = nrow(omega);
  
  om.old   = omega;
  df.old   = df
  psi.old  = df$psi

  om.new = omega;

  ## ---- ppsl
  
  mv.new = get.laplace.2(omega, y, XiL, df, DP, block.start, block.end, bsize)
  ## CL.new = t(chol(mv.new$V));
  U.new  = chol(mv.new$P);
  if (just.max) {
    ep.new = rep(0.0, nblocks*bsize);
  } else {
    ep.new = rnorm(nblocks*bsize);    
    ## om.new.sm = mv.new$m + CL.new %*% ep.new;
  }
  om.new.sm = mv.new$m + solve(U.new, ep.new);
  om.new.sm = as.numeric(om.new.sm);
  om.new[,inb.idc] = om.new.sm
  beta.new  = L %*% as.numeric(om.new);
  psi.new   = colSums(t(X) * as.numeric(beta.new));
  ## lpp.new   = -1 * sum(log(diag(CL.new))) - 0.5 * ep.new %*% ep.new;
  lpp.new   = sum(log(diag(U.new))) - 0.5 * ep.new %*% ep.new;

  df.new    = logit.llh(y, psi.new, n);

  ## cat("ep:", ep.new, "\n");
  ## cat("m:", mv.new$m, "\n");
  ## cat("d:", om.new.sm, "\n");
  
  ## cat("psi.new", psi.new, "\n");
  ## cat("psi.old", psi.old, "\n");
  ## cat("om.new", om.new, "\n");
  ## cat("beta.new", beta.new, "\n");
  
  ## ---- old

  mv.old = get.laplace.2(om.old, y, XiL, df.new, DP, block.start, block.end, bsize)
  ## CL.old = t(chol(mv.old$V));
  U.old     = chol(mv.old$P);
  om.old.sm = as.numeric(om.old[,inb.idc]);
  ## ep.old    = solve(CL.old, om.old.sm - mv.old$m);
  ep.old = U.old %*% (om.old.sm - mv.old$m);
  ## lpp.old   = -1 * sum(log(diag(CL.old))) - 0.5 * t(ep.old) %*% ep.old;
  lpp.old   = sum(log(diag(U.old))) - 0.5 * t(ep.old) %*% ep.old;

  llh.ratio    = sum(df.new$l0 - df.old$l0);

  om.new.num  = as.numeric(om.new);
  om.old.num  = as.numeric(om.old);
  
  lprior.ratio = -0.5 * t(om.new.num) %*% DP %*% om.new.num -
                 -0.5 * t(om.old.num) %*% DP %*% om.old.num
  
  lpp.ratio    = lpp.new - lpp.old

  ## cat("llh.ratio:", llh.ratio, "\n");
  ## cat("lprior.ratio", lprior.ratio, "\n");
  ## cat("lpp.ratio:", lpp.ratio, "\n");

  ## cat("llike.new:", sum(df.new$l0), "llike.old:", sum(df.old$l0), "diff:", llh.ratio, "\n")
  ## cat("lpp.new:", lpp.new, "lpp.old", lpp.old, "diff:", lpp.ratio, "\n");
  ## cat("lprior.ratio:", lprior.ratio, "\n");
  
  l.ratio   = llh.ratio + lprior.ratio - lpp.ratio
  
  omega  = om.old
  df     = df.old
  accept = 0;
  psi    = df.old$psi

  accept = log(runif(1)) < l.ratio
  accept = accept | just.max
  ## accept = TRUE

  if (accept) {
    omega  = om.new
    df     = df.new
    psi    = df.new$psi
  }
  
  out = list("df"=df, "l.ratio"=l.ratio, "psi"=psi, "omega"=omega, "accept"=accept);
}

draw.omega.2 <- function(y, X, omega, df, XiL, L, DP, starts, n=1, just.max=FALSE)
{
  ## Make sure to check starts.
  
  n.starts = length(starts);
  N = length(y);
  ends = c(starts[-1]-1, N);
  naccept = 0;
  ntotal = 0;

  l.ratio = rep(0, n.starts);

  for (i in 1:n.starts) {
    block.start = starts[i]
    block.end   = ends[i]
    ## cat(block.start, block.end, "\n");
    draw  = draw.omega.block.2(y, X, omega, df, XiL, L, DP, block.start, block.end, n, just.max)
    omega = draw$omega;
    df    = draw$df;

    l.ratio[i] = draw$l.ratio
    
    naccept = naccept + draw$accept
    ntotal = ntotal + 1
  }

  out = list("df"=df, "l.ratio"=l.ratio, "psi"=draw$psi,
    "naccept"=naccept, "ntotal"=ntotal, "omega"=omega);
  out
}

################################################################################

if (FALSE) {

  ## source("DynLogitOmegaBlockMH.R")
  N = 10
  P = 1

  phi = 0.9
  sig = 1.0
  n = 20
  C0 = diag(1,P)
  
  beta.t  = matrix(0, nrow=P, ncol=N);
  omega.t = matrix(rnorm(N*P, 0, sig), nrow=P, ncol=N);
  
  beta.1  = t(chol(C0)) %*% rnorm(P) + rnorm(P, 0, sig);
  beta.t[,1]  = beta.1;
  omega.t[,1] = beta.1;
  
  for (i in 2:N) {
    beta.t[,i] = phi * beta.t[,i-1] + omega.t[,i];
  }

  X = matrix(1, nrow=N, ncol=P)
  psi.t = colSums(t(X) * beta.t);
  p.t = exp(psi.t) / (1 + exp(psi.t))
  y = rbinom(N, n, prob=p.t);

  plot(y, ylim=c(0, n))
  lines(n*p.t, col=2);

  ## Make precision matrix
  D = matrix(0, N*P, N*P);
  A = D;
  for (i in 1:N) {
    didc = bidc.to.midc(i, P);
    D[didc, didc] = diag(sig,P);
    A[didc, didc] = diag(1, P);
    if (i == 1) {
      D[didc,didc] = diag(sig, P) + diag(phi,P) %*% C0 %*% diag(phi,P);
    } else {
      A[didc,didc-P] = -diag(phi,P);
    }
  }

  om.prec = solve(D);
  prec = t(A) %*% solve(D) %*% A;

  L = solve(A);
  
}

if (FALSE) {

  ## Assume X is design matrix.
  N  = nrow(X)
  P  = ncol(X)

  tDX = matrix(0, nrow=N, ncol=N*P);
  for (i in 1:N) {
    ii = (i-1)*P + 1;
    tDX[i,ii:(ii+P-1)] = X[i,]
  }

  XiL = tDX %*% L;
  
}

################################################################################

if (FALSE) {

  ## source("DynLogitOmegaBlockMH.R")
  df = logit.llh(y, rep(0, N), n);
  M = 4000
  block.start = 1;
  block.end   = N;
  ## block.end   = N;
  
  output1 <- list("omega"=array(0,dim=c(M, P, N)),
                 "psi"=matrix(0,nrow=M, ncol=N));
  ## beta = matrix(0, P, N);
  omega = omega.t
  omega[block.start:block.end] = 0;
  naccept = 0;

  beta = L %*% as.numeric(omega);
  psi  = colSums(t(X) * as.numeric(beta));
  df   = logit.llh(y, psi, n);

  ## source("DynLogitOmegaBlockMH.R")
  for (i in 1:M) {
    draw = draw.omega.block.2(y=y, X=X, omega=omega, df=df,
      XiL=XiL, L=L, DP=om.prec,
      block.start=block.start, block.end=block.end, n=n)
    omega = draw$omega
    df = draw$df;
    output1$omega[i,,] = omega
    output1$psi[i,]   = draw$psi;
    naccept = naccept +  draw$accept;
  }

  om.m1 = apply(output1$omega, c(2,3), mean);
  ## om.m1 = output1$omega[M,,]
  p.post.om1 = 1 / (1 + exp(-output1$psi));
  p.om1    = apply(p.post.om1, 2, mean)
  beta.om1 = matrix(L %*% as.numeric(om.m1), nrow=P, ncol=N);
  
  naccept
  beta.t
  beta.om1
 
  plot(y)
  lines(n*p.t, col=2);
  lines(n*p.om1, col=3);
  ## lines(n*p.m1, col=3, lty=2);
  
}

##------------------------------------------------------------------------------

if (FALSE) {

  ## source("DynLogitOmegaBlockMH.R")
  df = logit.llh(y, rep(0, N), n);
  M = 4000
  burn = 1000
  starts = c(1, 6)
  
  output2 <- list("omega"=array(0,dim=c(M, P, N)),
                 "psi"=matrix(0,nrow=M, ncol=N),
                 "l.ratio"=matrix(0,nrow=M,ncol=length(starts)));
  ## beta = matrix(0, P, N);
  omega = omega.t
  naccept = 0;

  naccept = 0;
  ntotal  = 0;

  set.seed(1000);
  
  ## source("DynLogitOmegaBlockMH.R")
  omega = matrix(0,P,N)
  beta = L %*% as.numeric(omega);
  psi  = colSums(t(X) * as.numeric(beta));
  df   = logit.llh(y, psi, n);
  start.time = proc.time();
  for (i in 1:(M+burn)) {
    draw = draw.omega.2(y, X, omega, df, XiL, L, DP=om.prec, starts, n=n, just.max=FALSE)
    omega = draw$omega
    df   = draw$df
    if (i > burn) {
      ii = i - burn
      output2$omega[ii,,] = omega;
      output2$psi[ii,]   = draw$psi;
      output2$l.ratio[ii,] = draw$l.ratio;
      naccept = naccept + draw$naccept
      ntotal  = ntotal  + draw$ntotal
    }
    if (i %% 1000 == 0) cat("Iter:", i, "arate:", naccept/ntotal, "\n");
  }
  total.time = proc.time() - start.time
  print(total.time);

  om.m2 = apply(output2$omega, c(2,3), mean);
  ## om.m2 = output2$omega[M,,]
  p.post.om2 = 1 / (1 + exp(-output2$psi));
  p.om2    = apply(p.post.om2, 2, mean)
  beta.om2 = matrix(L %*% as.numeric(om.m2), nrow=P, ncol=N);
  
  naccept
  beta.t
  beta.om2
 
  plot(y)
  lines(n*p.t, col=1);
  lines(n*p.om2, col=2, lty=2);
  lines(n*p.om1, col=3, lty=2);
  ## lines(n*p.m1, col=3, lty=3);
  
}

##------------------------------------------------------------------------------

if (FALSE) {
  
  ## source("DynLogitOmegaBlockMH.R"); source("ManualLoad.R")
  df = logit.llh(y, rep(0, N), n);
  M = 10
  burn = 0
  starts = c(1)
  
  output3 <- list("omega"=array(0,dim=c(M, P, N)),
                  "beta"=array(0,dim=c(M,P,N)),
                  "psi"=matrix(0,nrow=M, ncol=N),
                  "l.ratio"=matrix(0,nrow=M,ncol=length(starts)));
  ## beta = matrix(0, P, N);
  omega = omega.t
  naccept = 0;
  omega = matrix(0, P, N)
  
  tX = t(X);
  
  beta = L %*% as.numeric(omega);
  psi  = colSums(tX * as.numeric(beta));
  df   = logit.llh.3(y, psi, 0, n);
  ## df$pl2 = df$phi;
  beta = matrix(beta, 1, N);

  naccept = 0;
  ntotal  = 0;

  ## source("ManualLoad.R")
  ## draw  = draw.omega.mh(omega, beta, df, y, tX, ntrials=n, prior.prec=om.prec, phi, starts, just.max=TRUE)

  set.seed(1000);
  
  ## source("DynLogitOmegaBlockMH.R")
  reps = (M+burn)
  start.time = proc.time()
  for (i in 1:reps) {
    draw  = draw.omega.mh(omega, beta, df, y=y, tX=tX, ntrials=n, prior.prec=om.prec, phi, starts, just.max=TRUE)
    omega = draw$omega
    beta  = draw$beta
    df    = draw$llh
    if (i > burn) {
      ii = i - burn
      output3$omega[ii,,] = omega;
      output3$beta[ii,,] = beta;
      output3$psi[ii,]   = draw$llh$psi;
      naccept = naccept + draw$naccept
    }
    ## print(omega)
    if (i %% 1000 == 0) cat("Iter:", i, "naccept:", naccept, "\n");
  }
  total.time = proc.time() - start.time
  print(total.time);

  om.m3 = apply(output3$omega, c(2,3), mean);
  ## om.m3 = output3$omega[M,,]
  p.post.om3 = 1 / (1 + exp(-output3$psi));
  p.om3    = apply(p.post.om3, 2, mean)
  beta.om3 = matrix(L %*% as.numeric(om.m3), nrow=P, ncol=N);
  
  naccept
  beta.t
  beta.om3
 
  plot(y)
  lines(n*p.t, col=1);
  lines(n*p.om3, col=4);
  lines(n*p.om1, col=3, lty=2);
  ## lines(n*p.om2, col=3, lty=3);
  
}

################################################################################

if (FALSE) {

  ## source("DynLogitOmegaBlockMH.R")
  N = 10
  P = 1

  P.a = 2
  X.stc = matrix(rnorm(N*2), N, 2);
  ## alpha.t = c(-1, 1);
  alpha.t = c(0, 0);
  
  phi = 0.9
  sig = 1.0
  W   = sig^2;
  n = 20
  m0 = rep(0, P)
  C0 = diag(1,P)
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

  X = matrix(1, nrow=N, ncol=P)
  psi.dyn.t = colSums(t(X) * beta.t)
  psi.stc.t = X.stc %*% alpha.t
  psi.t = psi.dyn.t + psi.stc.t;
  p.t = exp(psi.t) / (1 + exp(psi.t))
  y = rbinom(N, n, prob=p.t);

  plot(y, ylim=c(0, n))
  lines(n*p.t, col=2);
  
}

if (FALSE) {

  ## starts = c(1,5)
  starts = c(1, 6);
  samp = 4000
  burn = 1000
  verbose = 1000
  just.max = FALSE

  set.seed(1000);
  
  source("DynLogitOmegaBlock.R")
  output4 <- dyn.logit.om(y=y, X.dyn=X, n=n, m0=m0, C0=C0,
                      samp=samp, burn=burn, verbose=verbose, starts=starts,
                      mu.m0=NULL, mu.P0=NULL,
                      phi.m0=NULL, phi.P0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      X.stc=NULL, m0.stc=m0.stc, C0.stc=C0.stc,
                      mu.true = 0, phi.true=phi, W.true=W,
                      alpha.true=NULL, beta.true = NULL, 
                      just.max=just.max)

  ## plot(beta.t[1,], type="l")
  ## lines(output4$beta[samp,1,], col=2)

  ## alpha.m4 = output4$alpha[samp,]
  ## beta.m4  = output4$beta[samp,,]
  ## psi.m4   = colSums(t(X) * beta.m4) + X.stc %*% alpha.m4
  psi.m4   = output4$psi
  p.m4     = exp(psi.m4) / (1 + exp(psi.m4))
  p.m4     = apply(p.m4, 2, mean)

  plot(y)
  lines(n * p.t, col=1)
  lines(n * p.m4, col=4, lty=4)

  ##------------------------------------------------------------------------------
  
  pg.C0 = diag(1, P + P.a); pg.C0[1:P.a,1:P.a] = C0.stc; pg.C0[P.a:(P.a+P),P.a:(P.a+P)] = C0
  pg.m0 = c(m0.stc, m0)

  pg.C0 = C0
  pg.m0 = m0

  source("DynLogitPG.R")
  out.pg <- dyn.logit.PG(y=y, X.dyn=X, n=rep(n, length(y)), X.stc=NULL,
                         samp=samp, burn=burn, verbose=verbose,
                         m.0=pg.m0, C.0=pg.C0,
                         mu.m0=NULL, mu.P0=NULL,
                         phi.m0=NULL, phi.P0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         beta.true=NULL, iota.true=NULL, w.true=NULL,
                         mu.true=0, phi.true=phi, W.true=W)
    
  ## beta.pg  = apply(out.pg$beta, c(2,3), mean)
  ## alpha.pg = apply(out.pg$alpha, 2, mean)
  ## psi.m5 = colSums(t(X) * beta.pg[,-1]) + X.stc %*% alpha.pg
  p.m5   = 1.0 / (1 + exp(-out.pg$psi))
  p.m5   = apply(p.m5, 2, mean)
  lines(n * p.m5, col=5, lty=5);

  ##------------------------------------------------------------------------------
  
 ## source("DynLogitOmegaBlockMH.R"); source("ManualLoad.R")
  df = logit.llh(y, rep(0, N), n);
  M = samp
  burn = burn
  
  output3 <- list("omega"=array(0,dim=c(M, P, N)),
                  "beta"=array(0,dim=c(M,P,N)),
                  "psi"=matrix(0,nrow=M, ncol=N),
                  "l.ratio"=matrix(0,nrow=M,ncol=length(starts)));
  ## beta = matrix(0, P, N);
  omega = omega.t
  naccept = 0;
  omega = matrix(0, P, N)
  
  tX = t(X);
  
  beta = L %*% as.numeric(omega);
  psi  = colSums(tX * as.numeric(beta));
  df   = logit.llh.3(y, psi, 0, n);
  ## df$pl2 = df$phi;
  beta = matrix(beta, 1, N);

  naccept = 0;
  ntotal  = 0;

  set.seed(1000);
  
  ## source("DynLogitOmegaBlockMH.R")
  reps = (M+burn)
  start.time = proc.time()
  for (i in 1:reps) {
    draw  = draw.omega.mh(omega, beta, df, y=y, tX=tX, ntrials=n, prior.prec=om.prec, phi, starts=starts, just.max=just.max)
    omega = draw$omega
    beta  = draw$beta
    df    = draw$llh
    if (i > burn) {
      ii = i - burn
      output3$omega[ii,,] = omega;
      output3$beta[ii,,] = beta;
      output3$psi[ii,]   = draw$llh$psi;
      naccept = naccept + draw$naccept
    }
    ## print(omega)
    if (i %% 1000 == 0) cat("Iter:", i, "naccept:", naccept, "\n");
  }
  total.time = proc.time() - start.time
  print(total.time);

  ## om.m3 = apply(output3$omega, c(2,3), mean);
  ## om.m3 = output3$omega[M,,]
  p.post.om3 = 1 / (1 + exp(-output3$psi));
  p.om3    = apply(p.post.om3, 2, mean)
  beta.om3 = matrix(L %*% as.numeric(om.m3), nrow=P, ncol=N);

  lines(n*p.om3, col=3, lty=3)
  
}
