## For benchmarking dynamic binomial logistic regression.

## library("BayesLogit")
library("coda")

source("Benchmark-Utilities.R");

source("DynLogitPG.R")
source("DynLogitdRUM.R")
source("DynLogitCUBS.R")
## source("CUBS-MH.R")

################################################################################
                                  ## SETUP ##
################################################################################

run <- list("tokyo"=FALSE,
            "low.2"=FALSE,
            "low.4"=FALSE,
            "high.2"=FALSE,
            "high.4"=FALSE,
            "allsynth"=FALSE)

write.dir = "./"

write.it = FALSE
plot.it  = FALSE
print.it = FALSE
read.it  = FALSE

methods = c("PG", "dRUM", "CUBS", "FS2010");

run.idc = 1:3

samp = 10000
burn  = 2000
verbose = 1000
ntrials = 2

################################################################################
                           ## Dyn Logit Benchmark ##
################################################################################

benchmark.dyn.logit <- function(y, X.dyn, n, X.stc=NULL, 
                                samp=1000, burn=100, ntrials=1, verbose=1000,
                                method="PG", var.names=c("beta"), dset.name="NULL",
                                m.0=NULL, C.0=NULL,
                                mu.m0=NULL, mu.P0=NULL,
                                phi.m0=NULL, phi.P0=NULL,
                                W.a0=NULL, W.b0=NULL,
                                beta.true=NULL, iota.true=NULL,
                                mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  T = length(y)
  ## n is total number of trials per observation - not vectorized, i.e. the same for each i.
  cat("Will benchmark", method[1], "using", dset.name, "dataset for variable(s)", var.names, "\n");
  
  sstat = list(); for (nm in var.names) sstat[[nm]] = list();
  arate    = rep(0, ntrials);
  ess.time = rep(0, ntrials);
  
  for(i in 1:ntrials) {

    if (method=="PG") {
      gb <- dyn.logit.PG(y=y, X.dyn=X.dyn, n=rep(n,T), X.stc=X.stc,
                         samp=samp, burn=burn, verbose=verbose,
                         m.0=m.0, C.0=C.0,
                         mu.m0=mu.m0, mu.P0=mu.P0,
                         phi.m0=phi.m0, phi.P0=phi.P0,
                         W.a0=W.a0, W.b0=W.b0,
                         beta.true=beta.true, iota.true=iota.true, w.true=NULL,
                         mu.true=mu.true, phi.true=phi.true, W.true=W.true)
      gb$a.rate = 1
    } else if (method=="dRUM") {
      gb <- dyn.logit.dRUM(y=y, X.dyn=X.dyn, n=n, X.stc=X.stc,
                         samp=samp, burn=burn, verbose=verbose,
                         m.0=m.0, C.0=C.0,
                         mu.m0=mu.m0, mu.P0=mu.P0,
                         phi.m0=phi.m0, phi.P0=phi.P0,
                         W.a0=W.a0, W.b0=W.b0,
                         z.true=NULL, r.true=NULL,
                         beta.true=beta.true, iota.true=iota.true,
                         mu.true=mu.true, phi.true=phi.true, W.true=W.true)
      gb$a.rate = 1
    } else if (method=="CUBS") {
      gb <- dyn.logit.CUBS(y=y, X.dyn=X.dyn, n=rep(n,T), m0=m.0, C0=C.0,
                           samp=samp, burn=burn, verbose=verbose,
                           mu.m0=mu.m0, mu.P0=mu.P0,
                           phi.m0=phi.m0, phi.P0=phi.P0,
                           W.a0=W.a0, W.b0=W.b0, X.stc=X.stc,
                           mu.true = mu.true, phi.true=phi.true, W.true=W.true)
      gb$a.rate = gb$ac.rate[samp]
    } else if (method=="FS2010") {
      gb <- dyn.logit.FS(y, X.dyn=X.dyn, n=n[1], X.stc=X.tc, ## assume n[i] = n[1].
                         samp=samp, burn=burn, verbose=verbose,
                         m.0=m.0, C.0=C.0,
                         mu.m0=mu.m0, mu.P0=mu.P0,
                         phi.m0=phi.m0, phi.P0=phi.P0,
                         W.a0=W.a0, W.b0=W.b0,
                         z.true=NULL, r.true=NULL,
                         beta.true=beta.true, iota.true=iota.true,
                         mu.true=mu.true, phi.true=phi.true, W.true=W.true)
    } else {
      print("Unknown method.")
      return(NA);
    }

    for (nm in var.names) {
      if (length(dim(gb[[nm]])) > 2) sstat[[nm]][[i]] = sum.stat.dyn(gb[[nm]], gb$ess.time[3], thin=1);
      if (length(dim(gb[[nm]])) ==2) sstat[[nm]][[i]] = sum.stat(gb[[nm]], gb$ess.time[3], thin=1);
    }

    arate[i]    = gb$a.rate
    ess.time[i] = gb$ess.time[3]

  }

  for (nm in var.names)
    sstat[[nm]] = simplify2array(sstat[[nm]]);
  
  out <- list("gb"=gb, "sstat"=sstat, "arate"=arate, "ess.time"=ess.time);

  out
}
################################################################################
                           ## BENCHMARK DATA SETS ##
################################################################################

##------------------------------------------------------------------------------
                              ## TOKYO RAINFALL ##
##------------------------------------------------------------------------------

if (run$tokyo) {

  ## source("Benchmark-DynLogit.R")
  tokyo = read.csv("Benchmark-DataSets/tokyo-rain.csv");
  tkrain = as.numeric(na.omit(as.numeric(as.matrix(tokyo))))

  T = length(tkrain)
  y = tkrain
  X = matrix(1, nrow=T, ncol=1)
  n = 2
  
  ## Prior
  b.m0 = -1.0;
  b.C0 = 1.0;
  W      = 0.1;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  bench.tokyo = list();
  if (read.it) load("bmark-tokyo.RData")
  
  ## source("Benchmark-DynLogit.R")
  for (i in run.idc) {
    ## source("Benchmark-DynLogit.R")
    nm = methods[i];
    bench.tokyo[[nm]] <- benchmark.dyn.logit(y, X.dyn=X, n=n, X.stc=NULL, 
                                     samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                     method=nm, var.names="beta", dset.name="Tokyo",
                                     m.0=b.m0, C.0=b.C0,
                                     W.a0=W.a0, W.b0=W.b0,
                                     mu.true=0.0, phi.true=1.0);
  }
  
  tokyo.table = setup.table.dyn(bench.tokyo, "beta")

  ## if (plot.it)  { plot.bench(pg, fs); plot.check.logit(y, X, n=n, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.tokyo, tokyo.table, file=file.path(write.dir, "bmark-tokyo.RData"))
  
}

if (FALSE) {

  ## par(mfrow=c(1,1))
  beta.mean = apply(plot.out$beta, c(2,3), mean);
  beta.95   = apply(plot.out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(plot.out$beta, c(2,3), function(x){quantile(x, 0.05)});

  ymin = min(beta.mean, plot.out$beta);
  ymax = max(beta.mean, plot.out$beta);

  plot(0:T, beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, beta.95[1,], col="pink")
  lines(0:T, beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  ## lines(0:T, beta[1,]);
  
  if (n[1] > 3) { points(1:T, log(y / (n-y)), cex=0.1) } else
  { points(1:T, (y-min(y)) / (max(y)-min(y)) * (ymax-ymin) + ymin, cex=0.1) }
  
  lines(0:T, plot.out$beta[100,1,], col=3);
  
}

if (FALSE) {

  par(mfrow=c(1,3))
  hist(plot.out$mu[,1])
  hist(plot.out$phi[,1])
  hist(plot.out$W[,1])

  apply(plot.out$mu,  2, mean)
  apply(plot.out$phi, 2, mean)
  apply(plot.out$W,   2, mean)
  
}

##------------------------------------------------------------------------------
                            ## P = 2, Low Corr X ##
##------------------------------------------------------------------------------

if (run$low.2)
{

  ## source("Benchmark-DynLogit.R")
  load("Benchmark-DataSets/DynLogit-synth-2-low-cor-X.RData")
  dset.name = "low-2";
  est.ar = "wout.ar"
  ## est.ar = "w.ar"
  filename = paste("bench-dynlogit-", dset.name, "-", est.ar, ".RData", sep="")
  
  T = length(y)
  X.dyn = X
  X.stc = matrix(1, nrow=T, ncol=1)
  P = ncol(X.dyn)
  
  ## Prior
  m0 = rep(0, P+1)
  C0 = diag(100, P+1)

  phi.m0 = rep(0.95, P);
  phi.V0 = rep(0.1,  P);
  W.guess = 0.1
  W.a0   = rep(300, P)
  W.b0   = W.a0 * W.guess
  mu.true = rep(0.0, P)
  phi.true = rep(phi.true, P)
  W.true = rep(W.true, P)
  
  bench.low.2 = list();
  if (read.it) load(filename)
  
  ## source("Benchmark-DynLogit.R")
  for (i in run.idc) {
    ## source("Benchmark-DynLogit.R")
    nm = methods[i];
    bench.low.2[[nm]] <- benchmark.dyn.logit(y, X.dyn=X.dyn, n=n[1], X.stc=X.stc, 
                                             samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                             method=nm, var.names="beta", dset.name=dset.name,
                                             m.0=m0, C.0=C0,
                                             phi.m0=phi.m0, phi.P0=1/phi.V0,
                                             W.a0=W.a0, W.b0=W.b0,
                                             mu.true=mu.true, phi.true=phi.true, W.true=W.true);
  }
  
  low.2.table = setup.table.dyn(bench.low.2, "beta")

  ## if (plot.it)  { plot.bench(pg, fs); plot.check.logit(y, X, n=n, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.low.2, low.2.table, file=file.path(write.dir, filename))
  
}

##------------------------------------------------------------------------------
                            ## P = 4, Low Corr X ##
##------------------------------------------------------------------------------

if (run$low.4)
{

  ## source("Benchmark-DynLogit.R")
  load("Benchmark-DataSets/DynLogit-synth-4-low-cor-X.RData")
  dset.name = "low-4";
  ## est.ar = "wout.ar"
  est.ar = "w.ar"
  filename = paste("bench-dynlogit-", dset.name, "-", est.ar, ".RData", sep="")
  
  T = length(y)
  X.dyn = X
  X.stc = matrix(1, nrow=T, ncol=1)
  P = ncol(X.dyn)
  
  ## Prior
  b.m0 = rep(0, P+1)
  b.C0 = diag(100, P+1)

  phi.m0 = rep(0.95, P);
  phi.V0 = rep(0.1,  P);
  W.guess = 0.1
  W.a0   = rep(300, P)
  W.b0   = W.a0 * W.guess
  mu.true = rep(0.0, P)
  phi.true = rep(phi.true, P)
  W.true = rep(W.true, P)
  
  bench.low.4 = list();
  if (read.it) load(file.name)
  
  ## source("Benchmark-DynLogit.R")
  for (i in run.idc) {
    ## source("Benchmark-DynLogit.R")
    nm = methods[i];
    bench.low.4[[nm]] <- benchmark.dyn.logit(y, X.dyn=X.dyn, n=n[1], X.stc=X.stc, 
                                             samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                             method=nm, var.names="beta", dset.name=dset.name,
                                             m.0=b.m0, C.0=b.C0,
                                             phi.m0=phi.m0, phi.P0=1/phi.V0,
                                             W.a0=W.a0, W.b0=W.b0,
                                             mu.true=mu.true, phi.true=NULL, W.true=NULL);
  }
  
  low.4.table = setup.table.dyn(bench.low.4, "beta")

  ## if (plot.it)  { plot.bench(pg, fs); plot.check.logit(y, X, n=n, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.low.4, low.4.table, file=file.path(write.dir, filename))
  
}

##------------------------------------------------------------------------------
                            ## P = 2, High Corr X ##
##------------------------------------------------------------------------------

if (run$high.2)
{

  ## source("Benchmark-DynLogit.R")
  load("Benchmark-DataSets/DynLogit-synth-2-high-cor-X.RData")
  dset.name = "high-2";
  est.ar = "wout.ar"
  ## est.ar = "w.ar"
  filename = paste("bench-dynlogit-", dset.name, "-", est.ar, ".RData", sep="")
  
  T = length(y)
  X.dyn = X
  X.stc = matrix(1, nrow=T, ncol=1)
  P = ncol(X.dyn)
  
  ## Prior
  b.m0 = rep(0, P+1)
  b.C0 = diag(100, P+1)

  phi.m0 = rep(0.95, P);
  phi.V0 = rep(0.1,  P);
  W.guess = 0.1
  W.a0   = rep(300, P)
  W.b0   = W.a0 * W.guess
  mu.true = rep(0.0, P)
  phi.true = rep(phi.true, P)
  W.true = rep(W.true, P)
  
  bench.high.2 = list();
  if (read.it) load(filename)
  
  ## source("Benchmark-DynLogit.R")
  for (i in run.idc) {
    ## source("Benchmark-DynLogit.R")
    nm = methods[i];
    bench.high.2[[nm]] <- benchmark.dyn.logit(y, X.dyn=X.dyn, n=n[1], X.stc=X.stc, 
                                             samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                             method=nm, var.names="beta", dset.name=dset.name,
                                             m.0=b.m0, C.0=b.C0,
                                             phi.m0=phi.m0, phi.P0=1/phi.V0,
                                             W.a0=W.a0, W.b0=W.b0,
                                             mu.true=mu.true, phi.true=phi.true, W.true=W.true);
  }
  
  high.2.table = setup.table.dyn(bench.high.2, "beta")

  ## if (plot.it)  { plot.bench(pg, fs); plot.check.logit(y, X, n=n, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.high.2, high.2.table, file=file.path(write.dir, filename))
  
}

##------------------------------------------------------------------------------
                            ## P = 4, High Corr X ##
##------------------------------------------------------------------------------

if (run$high.4)
{

  ## source("Benchmark-DynLogit.R")
  load("Benchmark-DataSets/DynLogit-synth-4-high-cor-X.RData")
  dset.name = "high-4";
  ## est.ar = "wout.ar"
  est.ar = "w.ar"
  filename = paste("bench-dynlogit-", dset.name, "-", est.ar, ".RData", sep="")
  
  T = length(y)
  X.dyn = X
  X.stc = matrix(1, nrow=T, ncol=1)
  P = ncol(X.dyn)
  
  ## Prior
  b.m0 = rep(0, P+1)
  b.C0 = diag(100, P+1)

  phi.m0 = rep(0.95, P);
  phi.V0 = rep(0.1,  P);
  W.guess = 0.1
  W.a0   = rep(300, P)
  W.b0   = W.a0 * W.guess
  mu.true = rep(0.0, P)
  phi.true = rep(phi.true, P)
  W.true = rep(W.true, P)
  
  bench.high.4 = list();
  if (read.it) load(filename)
  
  ## source("Benchmark-DynLogit.R")
  for (i in run.idc) {
    ## source("Benchmark-DynLogit.R")
    nm = methods[i];
    bench.high.4[[nm]] <- benchmark.dyn.logit(y, X.dyn=X.dyn, n=n[1], X.stc=X.stc, 
                                             samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                             method=nm, var.names="beta", dset.name=dset.name,
                                             m.0=b.m0, C.0=b.C0,
                                             phi.m0=phi.m0, phi.P0=1/phi.V0,
                                             W.a0=W.a0, W.b0=W.b0,
                                             mu.true=mu.true, phi.true=NULL, W.true=NULL);
  }
  
  high.4.table = setup.table.dyn(bench.high.4, "beta")

  ## if (plot.it)  { plot.bench(pg, fs); plot.check.logit(y, X, n=n, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.high.4, high.4.table, file=file.path(write.dir, filename))
  
}

##------------------------------------------------------------------------------
                                ## Synthetic ##
##------------------------------------------------------------------------------

if (run$allsynth)
{

  ## source("Benchmark-DynLogit.R")

  P = 2
  nt = 1
  corr.type = "low"
  est.ar = "with.ar"
  ## est.ar = "wout.ar"

  ## for (est.ar in c("wout.ar", "with.ar")) {
    ## for (P in c(2,4)) {
      for (nt in c(1,10)) {
        for (corr.type in c("low", "high")) {


  cat("Dyn Logit.  AR:", est.ar, "\n");
          
  dset.name = paste(corr.type, "-", P, "-n-", nt, sep="");
  source.file = paste("DynLogit-synth-", dset.name, ".RData", sep="")
  load(file.path("Benchmark-DataSets", source.file))

  filename = paste("bench-dynlogit-", dset.name, "-", est.ar, ".RData", sep="")
  
  T = length(y)
  X.dyn = X
  X.stc = matrix(1, nrow=T, ncol=1)
  P = ncol(X.dyn)
  
  ## Prior
  m0 = rep(0, P+1)
  C0 = diag(1, P+1)

  phi.m0 = rep(0.95, P);
  phi.V0 = rep(0.01,  P);
  
  W.guess = 0.1
  W.a0   = rep(10, P)
  W.b0   = W.a0 * W.guess
  mu.true = rep(0.0, P)

  phi.send = phi.true
  W.send   = W.true
  
  if (est.ar=="with.ar") {
    phi.send = NULL
    W.send = NULL
  }

  bench.synth = list();
  if (read.it) load(filename)
  
  ## source("Benchmark-DynLogit.R")
  for (i in run.idc) {
    ## source("Benchmark-DynLogit.R")
    nm = methods[i];
    bench.synth[[nm]] <- benchmark.dyn.logit(y, X.dyn=X.dyn, n=n[1], X.stc=X.stc, 
                                             samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                             method=nm,
                                             var.names=c("beta", "alpha", "phi", "W"),
                                             dset.name=dset.name,
                                             m.0=m0, C.0=C0,
                                             phi.m0=phi.m0, phi.P0=1/phi.V0,
                                             W.a0=W.a0, W.b0=W.b0,
                                             mu.true=mu.true, phi.true=phi.send, W.true=W.send);
  }
  
  synth.table = setup.table.dyn(bench.synth, "beta")

  ## if (plot.it)  { plot.bench(pg, fs); plot.check.logit(y, X, n=n, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.synth, synth.table, dset.name, file=file.path(write.dir, filename))

}}#}}
  
}


if (FALSE) {

  P = 2
  nt = 1
  corr.type = "low"
  est.ar = "with.ar"
  ## est.ar = "wout.ar"
  
  for (est.ar in c("wout.ar", "with.ar")) {
    ## for (P in c(2,4)) {
    for (nt in c(1,10)) {
      for (corr.type in c("low", "high")) {
        
        dset.name = paste(corr.type, "-", P, "-n-", nt, sep="");
        source.file = paste("DynLogit-synth-", dset.name, ".RData", sep="")
        bench.base  = paste("bench-dynlogit-", dset.name, "-", est.ar, sep="");
        bench.data  = paste(bench.base, ".RData", sep="")
        table.file  = paste("table.", bench.base, sep="");
      
        cat("AR:", est.ar, dset.name, "\n");        
        
        load(file.path("Benchmark-DataSets", source.file))
        load(file.path("Bench-Dyn-02", bench.data))
        
        the.table = synth.table$table
        the.table[3,2] = mean(bench.synth$CUBS$arate)
        
        write.table(the.table, file=table.file, row.names=TRUE, col.names=TRUE);
        
        par(mfrow=c(P+1,1))
        for (i in 1:P) {
          plot(beta[i,], col=1, type="l")
          ## plot(synth.table$ave.sstat[i,,1,1], type="l", col=2)
          lines(synth.table$ave.sstat[i,,1,1,drop=FALSE], col=2)
          lines(synth.table$ave.sstat[i,,1,2,drop=FALSE], col=3)
          ## lines(synth.table$ave.sstat[i,,1,3,drop=FALSE], col=4)
        }

        alpha.pg = mean(bench.synth$PG$gb$alpha);
        alpha.fs = mean(bench.synth$dRUM$gb$alpha);
        
        psi.pg = apply((synth.table$ave.sstat[,-1,1,1]) * t(X), 2, sum) + alpha.pg
        psi.fs = apply((synth.table$ave.sstat[,-1,1,2]) * t(X), 2, sum) + alpha.fs
        
        psi.max = max(psi)
        psi.min = min(psi)
        ymax = max(c(psi.pg, psi.max))
        ymin = min(c(psi.pg, psi.min))
        ytil = (y / nt) * (ymax - ymin) + ymin;
        
        plot(ytil, ylim=c(ymin,ymax))
        lines(psi, col=5)
        lines(psi.pg, col=2)
        lines(psi.fs, col=3)
        
        readline("<ENTER>")
        
      }}}
  
}
