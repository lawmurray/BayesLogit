

bidc.to.midc <- function(j, bsize, e=j) {
  bsize*(j-1) + 1:((e-j+1)*bsize);
}

logit.llh <- function(y, psi, n=1)
{
  epsi = exp(psi)
  p  = epsi / (1 + epsi);
  l0 = psi * y - n * log(1+epsi);
  l1 = y - n * p
  l2 = -n * (p - p^2)
  phi = psi * l2
  out = data.frame("l0"=l0, "l1"=l1, "l2"=l2, "phi"=phi, "psi"=psi, "pl2"=phi)
  out
}

logit.llh.2 <- function(y, X, beta, n=1)
{
  psi = colSums(t(X) * beta);
  epsi = exp(psi)
  p  = epsi / (1 + epsi);
  l0 = psi * y - n * log(1+epsi);
  l1 = y - n * p
  l2 = -n * (p - p^2)
  phi = psi * l2
  out = data.frame("l0"=l0, "l1"=l1, "l2"=l2, "phi"=phi, "psi"=psi, "pl2"=phi)
  out
}

logit.llh.3 <- function(y, psi.dyn, psi.stc, n=1)
{
  psi = psi.dyn + psi.stc
  epsi = exp(psi)
  p  = epsi / (1 + epsi);
  l0 = psi * y - n * log(1+epsi);
  l1 = y - n * p
  l2 = -n * (p - p^2)
  pl2 = psi * l2
  out = data.frame("psi"=psi, "psi.dyn"=psi.dyn, "psi.stc"=psi.stc, "l0"=l0, "l1"=l1, "l2"=l2, "pl2"=pl2)
  out
}

nb.mu.llh.3 <- function(y, psi.dyn, psi.stc, d=1)
{
  psi = psi.dyn + psi.stc
  epsi = exp(psi)
  bexpon = y + d
  p  = epsi / (d + epsi);
  l0 = psi * y - bexpon * log(d+epsi);
  l1 = y - bexpon * p
  l2 = -bexpon * (p - p^2)
  pl2 = psi * l2
  out = data.frame("psi"=psi, "psi.dyn"=psi.dyn, "psi.stc"=psi.stc, "l0"=l0, "l1"=l1, "l2"=l2, "pl2"=pl2)
  out
}

get.laplace.1 <- function(beta, Prec, X, df, block.start, block.end)
{
  ## Assume X is T x P in the form of a design matrix.
  ## Assume beta is P x T.
  ## Assume df is T x 5 df.
  bsize   = ncol(X)
  nblocks = block.end - block.start + 1;
  inb.idc = block.start:block.end
  inm.idc = (bsize * (block.start-1) + 1):(bsize * block.end);

  ## print(inm.idc);
  P0     = Prec[inm.idc, inm.idc,drop=FALSE]
  Om.11  = Prec[inm.idc, inm.idc,drop=FALSE]
  omb.12 = 0.0
  if (length(inm.idc) != nrow(Prec)) {
    Om.12  = Prec[inm.idc,-inm.idc,drop=FALSE]
    ## print(Om.12);
    beta.2 = as.numeric( beta[,-inm.idc] );
    ## print(beta.2)
    omb.12 = Om.12 %*% beta.2;
  }
  j = 0;
  for (i in inb.idc) {
    j = j + 1;
    diag.idx = 1:bsize + (j - 1)*bsize;
    x.i = as.numeric(X[i,])
    dxx.i = df$l2[i] * (x.i %*% t(x.i))
    #print(dxx.i);
    #print(Om.11[diag.idx,diag.idx])
    Om.11[diag.idx,diag.idx] = Om.11[diag.idx, diag.idx] - dxx.i
  }

  Gmat = X[inb.idc,] * df$l1 [inb.idc]
  Fmat = X[inb.idc,] * df$phi[inb.idc]
  G1 = as.numeric(t(Gmat))
  F1 = as.numeric(t(Fmat))
  ## cat("l1:", df$l1, "\n");
  ## cat("G1:", G1, "\n");
  ## cat("F1:", F1, "\n");
  b1 = G1 - F1

  V.post = solve(Om.11)
  m.post = - V.post %*% (omb.12 - b1);

  ##print(m.post)
  ##print(V.post)

  ## could also return a draw and log-dens calculation.
  ## or a precision and cholesky of precision.
  
  out = list("P0"=P0, "P"=Om.11, "V"=V.post, "m"=m.post);
  out
}

draw.beta.block.1 <- function(y, X, beta, df, prior.prec, block.start, block.end, n=1)
{
  ## Assume X is N x P like a design matrix.
  
  inb.idc = block.start:block.end;
  nblocks = block.end - block.start + 1;
  bsize   = nrow(beta);
  
  beta.old = beta;
  df.old   = df
  psi.old  = colSums(t(X) * beta.old);

  ppsl.beta = beta;

  ppsl.dist = get.laplace.1(beta.old, prior.prec, X, df.old, block.start, block.end);
  ppsl.L    = t(chol(ppsl.dist$V));
  ppsl.ep   = rnorm(nblocks*bsize)
  ppsl.beta.small = as.numeric(ppsl.dist$m + ppsl.L %*% ppsl.ep);
  ppsl.beta[,inb.idc] = ppsl.beta.small;
  ppsl.psi  = colSums(t(X) * ppsl.beta);
  ppsl.lpr  = -1 * sum(log(diag(ppsl.L))) - 0.5 * ppsl.ep %*% ppsl.ep;
  ## ppsl.lpr  = dmvnorm(x=ppsl.beta.small, mean=ppsl.dist$m, sigma=ppsl.dist$V, log=TRUE)
  
  ppsl.df   = logit.llh(y, ppsl.psi, n);

  dist.old  = get.laplace.1(ppsl.beta, prior.prec, X, ppsl.df, block.start, block.end);
  L.old     = t(chol(dist.old$V));

  beta.old.vec = as.numeric(beta.old[,inb.idc]);
  ep.old       = solve(L.old, beta.old.vec - dist.old$m);
  ## print(ep.old);
  lpr.old      = -1 * sum(log(diag(L.old))) - 0.5 * t(ep.old) %*% ep.old;
  ## lpr.old      = dmvnorm(beta.old.vec, dist.old$m, dist.old$V, log=TRUE)
  
  llh.ratio    = sum(ppsl.df$l0[inb.idc] - df.old$l0[inb.idc])
  ## lprior.ratio = -0.5 * t(ppsl.beta.small) %*% ppsl.dist$P0 %*% ppsl.beta.small -
  ##               -0.5 * t(beta.old.vec) %*% dist.old$P0 %*% beta.old.vec
  ppsl.beta.num = as.numeric(ppsl.beta);
  beta.num.old  = as.numeric(beta.old);
  lprior.ratio = -0.5 * t(ppsl.beta.num) %*% prior.prec %*% ppsl.beta.num -
                 -0.5 * t(beta.num.old) %*% prior.prec %*% beta.num.old
  lpr.ratio    = ppsl.lpr - lpr.old

  ## cat("llh.ratio:", llh.ratio, "\n");
  ## cat("lprior.ratio", lprior.ratio, "\n");
  ## cat("lpr.ratio:", lpr.ratio, "\n");

  l.ratio   = llh.ratio + lprior.ratio - lpr.ratio
  ## l.ratio = 1
  ## l.ratio = llh.ratio + lprior.ratio
  
  ## cat("l.ratio:", l.ratio, "\n");
  
  beta = beta.old
  df   = df.old
  accept = 0;
  psi  = psi.old
  
  if (log(runif(1)) < l.ratio) {
    accept = 1;
    beta   = ppsl.beta
    df     = ppsl.df
    psi    = ppsl.psi
  }
  
  out = list("df"=df, "l.ratio"=l.ratio, "psi"=psi, "beta"=beta, "accept"=accept);
}

draw.beta.1 <- function(y, X, beta, df, prior.prec, starts, n=1)
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
    draw = draw.beta.block.1(y, X, beta, df, prior.prec, block.start, block.end, n)
    beta = draw$beta;
    df   = draw$df;

    l.ratio[i] = draw$l.ratio
    
    naccept = naccept + draw$accept
    ntotal = ntotal + 1
  }

  out = list("df"=df, "l.ratio"=l.ratio, "psi"=draw$psi,
    "naccept"=naccept, "ntotal"=ntotal, "beta"=beta);
  out
}

################################################################################

if (FALSE) {

  ## source("DynLogitBetaBlockMH.R")
  N = 500
  P = 1

  phi = 0.9
  sig = 1.0
  n = 20
  C0 = diag(1,P)
  
  beta.t = matrix(0, nrow=P, ncol=N);
  for (i in 2:N) {
    beta.t[,i] = phi * beta.t[,i-1] + rnorm(P, 0, sig);
  }

  X = matrix(1, nrow=N, ncol=P)
  psi.t = colSums(t(X) * beta.t);
  p.t = exp(psi.t) / (1 + exp(psi.t))
  y = rbinom(N, n, prob=p.t);

  plot(y)
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

  prec = t(A) %*% solve(D) %*% A;

}

if (FALSE) {

  ## source("DynLogitBetaBlockMH.R")
  df = logit.llh(y, rep(0, N), n);
  M = 1000
  block.start = 1;
  block.end   = N;
  ## block.end   = N;
  
  output1 <- list("beta"=array(0,dim=c(M, P, N)),
                 "psi"=matrix(0,nrow=M, ncol=N));
  ## beta = matrix(0, P, N);
  beta = beta.t
  beta[block.start:block.end] = 0;
  naccept = 0;
    
  for (i in 1:M) {
    if (i %% 10 == 0) { cat("Iter:", i, "\n"); }
    draw = draw.beta.block.1(y, X, beta, df, prec, block.start, block.end, n)
    beta = draw$beta
    df = draw$df;
    output1$beta[i,,] = beta;
    output1$psi[i,]   = draw$psi;
    naccept = naccept +  draw$accept;
  }

  ## beta.m1 = apply(output1$beta, c(2,3), mean);
  beta.m1 = output$beta[M,,]
  p.post1 = 1 / (1 + exp(-output1$psi));
  p.m1    = apply(p.post1, 2, mean)
  
  naccept
  beta.t
  beta.m1
 
  plot(y)
  lines(n*p.t, col=2);
  lines(n*p.m1, col=3);
  lines(n*p.pg.m, col=4);
  
}


## Manual 1
if (FALSE) {

  ## source("DynLogitBetaBlockMH.R");
  M = 30
  inb.idc = block.start:block.end;
  beta = beta.t
  beta[inb.idc] = 1.0;
  beta.newton = matrix(0, nrow=M, ncol=N);
  for (i in 1:M) {
    df = logit.llh.2(y, X, beta, n);
    Om = prec - diag(df$l2, N)
    b  = df$l1 - df$psi * df$l2
    Om.11 = Om[inb.idc, inb.idc]
    Om.12 = Om[inb.idc, -inb.idc]
    V.11 = solve(Om.11)
    beta[inb.idc] = -V.11 %*% (Om.12 %*% beta[-inb.idc] - b[inb.idc]);
    beta.newton[i,] = beta;
  }

  beta.mode = beta
  psi.mode = colSums(t(X) * beta.mode);
  p.mode = 1 / (1 + exp(-psi.mode));
  
  plot(y)
  lines(n*p.t, col=2);
  lines(n*p.mode, col=3, lty=2);
  
}

## Manual 2
if (FALSE) {

  ## source("DynLogitBetaBlockMH.R");
  M = 10
  inb.idc = block.start:block.end;
  beta = beta.t
  beta[inb.idc] = 1.0;

  beta = rep(0, 10)
  beta.newton = matrix(0, nrow=M, ncol=N);
  for (i in 1:M) {
    print(beta)
    df = logit.llh.2(y, X, beta, n);
    Om = prec - diag(df$l2, N)
    b  = df$l1 - df$phi
    V.11 = solve(Om)
    beta = as.numeric(V.11 %*% b);
    beta.newton[i,] = beta;
  }

  print(beta)

  beta.mode = beta
  psi.mode = colSums(t(X) * beta.mode);
  p.mode = 1 / (1 + exp(-psi.mode));
  
  plot(y)
  lines(n*p.t, col=2);
  lines(n*p.mode, col=3, lty=2);

  ## ------------

  logit.llh <- function(beta, y, X, n, prec) {
    beta = as.numeric(beta)
    psi  = colSums(t(X) * beta);
    l    = - 0.5 * t(beta) %*% prec %*% beta + sum(psi*y - n * log(1+exp(psi)))
    l
  }

  logit.grad <- function(beta, y, X, n, prec) {
    beta = as.numeric(beta)
    psi  = colSums(t(X) * beta);
    g    = - 0.5 * prec %*% beta + n * exp(psi) / (1+exp(psi))
    g
  }

  beta = rep(0,N)
  oout = optim(beta, fn=logit.llh, y=y, X=X, n=n, prec=prec, method="Nelder-Mead",
    control=list("fnscale"=-1, "maxit"=10000));

  beta.mode2 = oout$par
  psi.mode2 = colSums(t(X) * beta.mode2);
  p.mode2 = 1 / (1 + exp(-psi.mode2));

  lines(n*p.mode2, col=5)

  ## -------------

  ## source("Metropolis.R")
  Xblah = diag(as.numeric(X), N)
  out.mh = ind.metropolis(y, Xblah, rep(n,length(y)),
    m0=matrix(0,N,1), P0=array(prec,dim=c(N,N,1)),
    m=beta.mode, V=V.11,
    samp=10000, burn=100, tune = 1.0, verbose=1000, df=Inf)

  beta.m3 = colMeans(out.mh$beta)
  psi.m3 = colSums(t(X) * beta.m3);
  p.mode3 = 1 / (1 + exp(-psi.m3));

  lines(n*p.mode3, col=6)
  
  
}

if (FALSE) {

  ## source("DynLogitBetaBlockMH.R")
  df = logit.llh(y, rep(0, N), n);
  M = 1000
  burn = 0
  ## starts = c(1,2);
  ## starts = 1:9
  starts = seq(1,491, 10);

  output <- list("beta"=array(0,dim=c(M, P, N)),
                 "psi"=matrix(0,nrow=M, ncol=N),
                 "l.ratio"=matrix(0,nrow=M,ncol=length(starts)));
  naccept = 0;
  ntotal  = 0;
  beta = matrix(0, nrow=P, ncol=N)
  psi = colSums(t(X) * beta)
  df = logit.llh(y, psi, n)
  
  for (i in 1:(M+burn)) {
    if (i %% 10 == 0) { cat("Iter:", i, "\n"); }
    draw = draw.beta.1(y, X, beta, df, prec, starts, n)
    beta = draw$beta
    df   = draw$df
    if (i > burn) {
      ii = i - burn
      output$beta[ii,,] = beta;
      output$psi[ii,]   = draw$psi;
      output$l.ratio[ii,] = draw$l.ratio;
      naccept = naccept + draw$naccept
      ntotal  = ntotal  + draw$ntotal
    }
  }

  beta.m = apply(output$beta, c(2,3), mean);
  p.post = 1 / (1 + exp(-output$psi));
  p.m    = apply(p.post, 2, mean)
  
  naccept
  ntotal
  #beta.t
  #beta.m
 
  plot(y)
  lines(n*p.t, col=2);
  lines(n*p.m, col=3);
  lines(n*p.mode, col=5, lty=4);
  lines(n*p.pg.m, col=4, lty=2);
  
}

if (FALSE ){

  M = 5000
  ## source("DynLogitPG.R")
  out.pg <- dyn.logit.PG(y=y, X.dyn=X, n=rep(n, length(y)), X.stc=NULL,
                         samp=M, burn=100, verbose=1000,
                         m.0=rep(0,P), C.0=C0,
                         beta.true=NULL, iota.true=0.0, w.true=NULL,
                         mu.true=0.0, phi.true=rep(phi,P), W.true=rep(sig^2,P))

  out.pg$psi = matrix(0, M, N);
  for (i in 1:M) {
    out.pg$psi[i,] = colSums(t(X) * out.pg$beta[i,,-1])
  }
  
  beta.pg.m = apply(out.pg$beta, c(2,3), mean);
  p.pg.post = 1 / (1 + exp(-out.pg$psi));
  p.pg.m    = apply(p.pg.post, 2, mean)
  

  beta.t
  beta.pg.m
 
  # plot(y)
  # lines(n*p, col=2);
  lines(n*p.pg.m, col=4);
  
  
}
