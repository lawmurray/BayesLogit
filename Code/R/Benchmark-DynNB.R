## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

## library("BayesLogit")
library("coda")

source("Benchmark-Utilities.R")

source("DynNBPG.R")
source("DynNBFS-2009.R")
source("DynNBCUBS.R")
source("DynNBOmegaBlock.R")

################################################################################
                                  ## SETUP ##
################################################################################

run <- list("flu"   =FALSE,
            "synth1"=FALSE,
            "synth2"=FALSE,
            "synth3"=FALSE,
            "allsynth"=FALSE,
            "allview"=FALSE)

write.dir = "Bench-Dyn-06" # Used in file.path

write.it = TRUE
plot.it  = FALSE
read.it = FALSE

methods = c("PG", "FS", "CUBS", "OmegaBlock")

## The methods to run from the sequence above.
run.idc = 1:3

samp = 10000
burn = 2000
verbose = 1000
ntrials = 10 ## Ambiguous languge.  Repetitions of MCMC.

options = list("just.max"=FALSE, "starts"=1);

################################################################################
                             ## Dyn NB Benchmark ##
################################################################################

benchmark.dyn.NB <- function(y, X.dyn, X.stc=NULL,
                             samp=1000, burn=100, ntrials=1, verbose=100,
                             method="PG", var.names=c("beta"), dset.name="NULL",
                             m.0=NULL, C.0=NULL,
                             mu.m0=NULL, mu.P0=NULL,
                             phi.m0=NULL, phi.P0=NULL,
                             W.a0=NULL, W.b0=NULL,
                             d.true=NULL,
                             beta.true=NULL, iota.true=NULL,
                             mu.true=NULL, phi.true=NULL, W.true=NULL, options=NULL)
{
  ## n is total number of trials per observation - not vectorized, i.e. the same for each i.
  cat("Will benchmark", method[1], "using", dset.name, "dataset for variable(s)", var.names, "\n");

  sstat = list(); for (nm in var.names) sstat[[nm]] = list();
  arate    = rep(0, ntrials);
  ess.time = rep(0, ntrials);
  gb = list();
  i  = 0
  num.errors = 0
  
  while (i < ntrials) {
    i = i + 1
    
    if (method=="PG") {
      gb <- dyn.NB.PG(y=y, X.dyn=X.dyn, X.stc=X.stc,
                      samp=samp, burn=burn, verbose=verbose,
                      m.0=m.0, C.0=C.0,
                      mu.m0=mu.m0, mu.P0=mu.P0,
                      phi.m0=phi.m0, phi.P0=phi.P0,
                        W.a0=W.a0, W.b0=W.b0,
                      d.true=d.true, w.true=NULL,
                      beta.true=beta.true, iota.true=iota.true,
                      mu.true=mu.true, phi.true=phi.true, W.true=W.true);
      gb$beta = gb$beta[,,-1,drop=FALSE]
      gb$a.rate = 1
    } else if (method=="FS") {
      gb <- dyn.NB.FS(y=y, X.dyn=X.dyn, X.stc=X.stc,
                      samp=samp, burn=burn, verbose=verbose,
                      m.0=m.0, C.0=C.0,
                      mu.m0=mu.m0, mu.P0=mu.P0,
                      phi.m0=phi.m0, phi.P0=phi.P0,
                      W.a0=W.a0, W.b0=W.b0,
                      d.true=d.true, lambda.true=NULL, r.true=NULL,
                      beta.true=beta.true, iota.true=iota.true,
                      mu.true=mu.true, phi.true=phi.true, W.true=W.true)
      gb$beta = gb$beta[,,-1,drop=FALSE]
      gb$a.rate=1
    } else if (method=="CUBS") {
      gb <- dyn.NB.CUBS(y=y, X.dyn=X.dyn, m0=m.0, C0=C.0,
                        samp=samp, burn=burn, verbose=verbose,
                        mu.m0=mu.m0, mu.P0=mu.P0,
                        phi.m0=phi.m0, phi.P0=phi.P0,
                        W.a0=W.a0, W.b0=W.b0, X.stc=X.stc,
                        mu.true = mu.true, phi.true=phi.true, W.true=W.true, d.true=d.true)
      gb$beta = gb$beta[,,-1,drop=FALSE]
      gb$a.rate = gb$ac.rate[samp]
    } else if (method=="OmegaBlock") {
      m0.stc = NULL
      C0.stc = NULL
      m0.dyn = m.0
      C0.dyn = C.0
      if (!is.null(X.stc)) {
        P.a = ncol(X.stc)
        P.b = ncol(X.dyn)
        idc.a = 1:P.a
        idc.b = (P.a+1):(P.a+P.b)
        m0.stc = m.0[idc.a]
        C0.stc = C.0[idc.a, idc.a]
        m0.dyn = m.0[idc.b]
        C0.dyn = C.0[idc.b,idc.b]
      }
      gb <- dyn.nb.om(y=y, X.dyn=X.dyn, m0=m0.dyn, C0=C0.dyn,
                      samp=samp, burn=burn, verbose=verbose, starts=options$starts,
                      mu.m0=mu.m0, mu.P0=mu.P0,
                      phi.m0=phi.m0, phi.P0=phi.P0,
                      W.a0=W.a0, W.b0=W.b0,
                      X.stc=X.stc, m0.stc=m0.stc, C0.stc=C0.stc,
                      mu.true = mu.true, phi.true=phi.true, W.true=W.true,
                      alpha.true=NULL, beta.true = NULL, d.true=d.true,
                      just.max=options$just.max)
      gb$options = options
    } else {
      print("Unknown method.")
      return(NA);
    }

    if (gb$error) {
      num.errors = num.errors + 1
      i = i-1
      cat("Error.  Dump:\n", gb$dump, "\n");
    }
    else {
    
      for (nm in var.names) {
        if (length(dim(gb[[nm]])) > 2) sstat[[nm]][[i]] = sum.stat.dyn(gb[[nm]], gb$ess.time[3], thin=1);
        if (length(dim(gb[[nm]])) ==2) sstat[[nm]][[i]] = sum.stat(gb[[nm]], gb$ess.time[3], thin=1);
      }
      
      arate[i]    = gb$a.rate[1]
      ess.time[i] = gb$ess.time[3]
      
    }

    if (num.errors > 10) return(NA)
  }

  for (nm in var.names)
    sstat[[nm]] = simplify2array(sstat[[nm]]);
  
  out <- list("gb"=gb, "sstat"=sstat, "arate"=arate, "ess.time"=ess.time);

  out
} ## benchmark.dyn.NB

################################################################################
                           ## BENCHMARK DATA SETS ##
################################################################################

##------------------------------------------------------------------------------
                                 ## FLU DATA ##
##------------------------------------------------------------------------------

if (run$flu) {

  flu = read.csv("../../fludata/flu5years.csv", header=TRUE)

  ## Can use any of these
  y1 = flu$"A.H3."
  y2 = flu$"A.H1."
  y3 = flu$B
  y4 = flu$Total.number.positive
  
  ## plot(y1)
  
  rawX = flu[,10:109]
  N = nrow(rawX)
  T = length(y1)
  
  ## Naive PCA
  ## Create a reduced-dimension representation of the search terms
  NPC = 4 ## Can play with this
  mypc = princomp(rawX)
  X = mypc$scores[,1:NPC]
  
  ## plot(X[,2])
  
  A = rawX - rep(1, N) %*% t(colMeans(rawX));
  sv = svd(A)
  Z = sv$u %*% diag(sv$d);
  cor(mypc$scores, Z)[1:4,1:4];
  ## you won't lose rank unless a column vector = alpha 1.
  
  ## So if we want to do a dynamic regression with intercept then we should use
  ## X.dyn = as.matrix(Z[,1:NPC])
  ## X.stc = matrix(1, nrow=T, ncol=1)
  ## X.stc = NULL
  ## Evolving intercept
  X.dyn = matrix(1, nrow=T, ncol=1)
  X.stc = Z[,1:NPC];

  N.b = ncol(X.dyn)
  if (!is.null(X.stc)) N.a = ncol(X.stc) else N.a = 0
  N   = N.b + N.a
  
  ## Prior
  mu.m0  = rep(0.0, N.b)
  mu.P0  = rep(0.1, N.b)
  phi.m0 = rep(0.95, N.b)
  phi.P0 = rep(100, N.b)
  m0   = rep(0.0, N);
  C0   = diag(3.0, N);
  W.gs   = 0.01
  W.a0   = rep(100, N.b);
  W.b0   = W.a0 * W.gs;

  bench.flu = list();
  
  ## source("Benchmark-DynNB.R")
  for (i in 1:3) {
    ## source("Benchmark-DynNB.R")
    nm = methods[i];
    bench.flu[[nm]] <- benchmark.dyn.NB(y=y1, X.dyn=X.dyn, X.stc=X.stc, 
                                        samp=samp, burn=burn, ntrials=1, verbose=verbose,
                                        method=nm, var.names="beta", dset.name="Flu",
                                        m.0=m0, C.0=C0,
                                        mu.m0=mu.m0, mu.P0=mu.P0,
                                        phi.m0=phi.m0, phi.P0=phi.P0,
                                        W.a0=W.a0, W.b0=W.b0, d.true=NULL);
    ## mu.true=rep(0.0, N.b), phi.true=rep(1, N.b)
  }
  
  flu.table = setup.table.dyn(bench.flu, "beta")
  
}

##------------------------------------------------------------------------------
                               ## SYNTHETIC 1 ##
##------------------------------------------------------------------------------

if (run$synth1) {

  load("Benchmark-DataSets/DynNB-synth-1.RData")

  y = dyn.nb.1$y;
  X.dyn = dyn.nb.1$X;
  P = 1
  N.y = length(y)
  
  options$starts=seq(1, N.y, 2)

  ## Prior
  b.m0 = 3.0;
  b.C0 = 3.0;
  W    = 0.1;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  bench.synth1 = list();
  
  ## source("Benchmark-DynNB.R")
  for (i in run.idc) {
    ## source("Benchmark-DynNB.R")
    nm = methods[i];
    bench.synth1[[nm]] <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=NULL, 
                                     samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                     method=nm, var.names="beta", dset.name="Synth1",
                                     m.0=b.m0, C.0=b.C0,
                                     W.a0=W.a0, W.b0=W.b0,
                                     mu.true=0.0, phi.true=1.0, d.true=4, options=options);
  }
  
  synth1.table = setup.table.dyn(bench.synth1, "beta")
 
  if (plot.it)  { plot.bench(pg, fs); plot.check.NB(y, X.dyn, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.synth1, synth1.table, file=file.path(write.dir, "bmark-synth1.RData"))
  
}

##------------------------------------------------------------------------------
                                 ## SYNTH 2 ##
##------------------------------------------------------------------------------

if (run$synth2) {

  load("Benchmark-DataSets/DynNB-synth-2.RData")

  y = dyn.nb.2$y;
  X.dyn = dyn.nb.2$X;

  ## Prior
  b.m0 = 3.0;
  b.C0 = 3.0;
  W    = 0.1;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  bench.synth2 = list();
  
  ## source("Benchmark-DynNB.R")
  for (i in run.idc) {
    ## source("Benchmark-DynNB.R")
    nm = methods[i];
    bench.synth2[[nm]] <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=NULL, 
                                     samp=samp, burn=burn, ntrials=1, verbose=verbose,
                                     method=nm, var.names="beta", dset.name="Synth2",
                                     m.0=b.m0, C.0=b.C0,
                                     W.a0=W.a0, W.b0=W.b0,
                                     mu.true=0.0, phi.true=1.0, options=options);
  }
  
  synth2.table = setup.table.dyn(bench.synth2, "beta")
 
  if (plot.it)  { plot.bench(pg, fs); plot.check.NB(y, X.dyn, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.synth2, synth2.table, file=file.path(write.dir, "bmark-synth2.RData"))

}

##------------------------------------------------------------------------------
                                 ## SYNTH 3 ##
##------------------------------------------------------------------------------

if (run$synth3) {

  load("Benchmark-DataSets/nb.synth3.RData")

  y = y;
  X.dyn = X.dyn;
  Ny = length(y)

  options$starts = c(1, Ny, 50)
  
  ## Prior
  b.m0 = 0.0;
  b.C0 = 3.0;
  W    = 0.1;
  W.a0   = 10;
  W.b0   = W.a0 * W;
  mu.m0 = 0.0
  mu.P0 = 1.0
  phi.m0 = 0.95
  phi.P0 = 100

  bench.synth3 = list();
  
  ## source("Benchmark-DynNB.R")
  for (i in run.idc) {
    ## source("Benchmark-DynNB.R")
    nm = methods[i];
    bench.synth3[[nm]] <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=NULL, 
                                     samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                     method=nm, var.names="beta", dset.name="Synth3",
                                     m.0=b.m0, C.0=b.C0,
                                     W.a0=W.a0, W.b0=W.b0,
                                     mu.true=0.0, phi.true=1.0, d.true = d, options=options);
  }
  
  synth3.table = setup.table.dyn(bench.synth3, "beta")
 
  if (write.it) save(bench.synth3, synth3.table, file=file.path(write.dir, "bmark-synth3.RData"))

}

##------------------------------------------------------------------------------
                              ## GENERIC SYNTH ##
##------------------------------------------------------------------------------

if (run$allsynth)
{

  ## source("Benchmark-DynNB.R")

  sstats = list()
  tables = list()
  ids    = list()
  iter   = 0
  
  P = 2
  nb.mean = 100
  corr.type = "low"
  ## est.ar = "with.ar"
  est.ar = "wout.ar"
  
  for (est.ar in c("wout.ar", "with.ar")) {
    ## for (P in c(2,4)) {
      for (nb.mean in c(10,100)) {
        for (corr.type in c("low", "high")) {

  iter = iter + 1
  cat("AR:", est.ar, "\n");
          
  dset.name = paste(corr.type, "-", P, "-mu-", nb.mean, sep="");
  source.file = paste("DynNB-synth-", dset.name, ".RData", sep="")
  load(file.path("Benchmark-DataSets", source.file))

  filename = paste("bench-dynnb-", dset.name, "-", est.ar, ".RData", sep="")
  
  T = length(y)
  X.dyn = X
  X.stc = matrix(1, nrow=T, ncol=1)
  P = ncol(X.dyn)

  options$starts = seq(1,T,50)
  
  ## Prior
  m0 = rep(0, P+1)
  C0 = diag(100, P+1)

  phi.m0 = rep(0.95, P);
  phi.V0 = rep(0.1,  P);
  W.guess = 0.1
  W.a0   = rep(300, P)
  W.b0   = W.a0 * W.guess
  mu.true = rep(0.0, P)

  if (est.ar=="with.ar") {
    phi.true = NULL
    W.true = NULL
  }    
  
  bench.synth = list();
  if (read.it) load(file.path(write.dir, filename))
  
  ## source("Benchmark-DynNB.R")
  for (i in run.idc) {
    ## source("Benchmark-DynNB.R")
    nm = methods[i];
    bench.synth[[nm]] <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=X.stc, 
                                          samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                          method=nm, var.names=c("beta", "alpha", "phi", "W"),
                                          dset.name=dset.name,
                                          m.0=m0, C.0=C0,
                                          phi.m0=phi.m0, phi.P0=1/phi.V0,
                                          W.a0=W.a0, W.b0=W.b0,
                                          mu.true=mu.true, phi.true=phi.true, W.true=W.true, d.true=d.true, options=options);
  }
  
  synth.table = setup.table.dyn(bench.synth, "beta")  
  
  ## if (plot.it)  { plot.bench(pg, fs); plot.check.logit(y, X, n=n, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.synth, synth.table, dset.name, file=file.path(write.dir, filename))

  sstats[[iter]] = synth.table$ave.sstat
  tables[[iter]] = synth.table$table
  ids[[iter]]    = filename
  
}}}#}
  
}

## Write tables.
if (FALSE) {

  for (i in 1:length(ids)) {
    tempname = sub("RData", "table", ids[[i]])
    write.table(tables[[i]], tempname)
  }
  
}

if (run$allview) {

  P = 2
  nb.mean = 10
  corr.type = "high"
  ## est.ar = "with.ar"
  est.ar = "wout.ar"
  
  for (est.ar in c("wout.ar", "with.ar")) {
    ## for (P in c(2,4)) {
    for (nb.mean in c(10,100)) {
      for (corr.type in c("low", "high")) {
        
        cat("AR:", est.ar, "P:", P, "nb.mean:", nb.mean, "corr:", corr.type, "\n");
        
        dset.name = paste(corr.type, "-", P, "-mu-", nb.mean, sep="");
        source.file = paste("DynNB-synth-", dset.name, ".RData", sep="")
        bench.base  = paste("bench-dynnb-", dset.name, "-", est.ar, sep="");
        bench.data  = paste(bench.base, ".RData", sep="")
        table.file  = paste("table.", bench.base, sep="");
                  
        load(file.path("Benchmark-DataSets", source.file))
        load(file.path(write.dir, bench.data))

        ## cat("phi.true", phi.true, "\n");
        ## cat("W.true:", W.true, "\n");
        
        the.table = synth.table$table
        the.table[3,2] = mean(bench.synth$CUBS$arate) ## CUBS
        
        write.table(the.table, file=table.file, row.names=TRUE, col.names=TRUE);

        par(mfrow=c(P+1,1))
        for (i in 1:P) {
          plot(beta[i,], col=1, type="l")
          ## plot(synth.table$ave.sstat[i,,1,1], type="l", col=2)
          lines(synth.table$ave.sstat[i,,1,1], col=2, lty=2)
          lines(synth.table$ave.sstat[i,,1,2], col=3, lty=3)
          lines(synth.table$ave.sstat[i,,1,3], col=4, lty=4) ## CUBS

        }

        alpha.pg = mean(bench.synth$PG$gb$alpha);
        alpha.fs = mean(bench.synth$FS$gb$alpha);
        alpha.cb = mean(bench.synth$CUBS$gb$alpha); ## CUBS

        lmean.pg = apply((synth.table$ave.sstat[,,1,1]) * t(X), 2, sum) + alpha.pg
        lmean.fs = apply((synth.table$ave.sstat[,,1,2]) * t(X), 2, sum) + alpha.fs
        lmean.cb = apply((synth.table$ave.sstat[,,1,3]) * t(X), 2, sum) + alpha.cb ## CUBS
        
        plot(y, cex=0.5)
        lines(exp(log.mean))
        lines(exp(lmean.pg), col=2, lty=2)
        lines(exp(lmean.fs), col=3, lty=3)
        lines(exp(lmean.cb), col=4, lty=4) ## CUBS

        readline("<ENTER>")
      }}}
  
}

################################################################################
#-------------------------------------------------------------------------------
                            ## GENERATE AND TEST ##
#-------------------------------------------------------------------------------
################################################################################

if (FALSE)
{

  T = 500
  P = 2
  corr.type = "low"
  nb.mean = 24

  ## for (P in c(2,4)) {
  ##   for (corr.type in c("low", "high")) {
  ##     for (nb.mean in c(10, 100)) {
  
  c  = 0.5
  d.true  = 4
  marg.V   = 5 / sqrt(P) * c
  phi.true = rep(0.95, P)

  W.true = marg.V * (1 - phi.true^2)

  beta = matrix(0, nrow=P, ncol=T+1)
  beta[,1] = 0
  for (i in 2:(T+1))
    beta[,i] = phi.true * (beta[,i-1]) + rnorm(P, 0, sqrt(W.true))
  
  xgrid = seq(-1, 1, length.out=T)
  
  tX = matrix(0, nrow=P, ncol=T);
  if (corr.type=="low")  freq = c(1, 2, 3, 4)
  if (corr.type=="high") freq = c(1, 1.1, 1.2, 1.3)
  for (i in 1:P)
    tX[i,] = cos(freq[i] * pi * xgrid);

  tX = tX / sqrt(P) * (1-c)
  X  = t(tX)
  
  log.mean = log(nb.mean) + colSums(beta[,-1] * tX)
  psi = log.mean - log(d.true)
  p.success  = 1 / (1 + exp(-psi))
  y   = rnbinom(T, d.true, 1-p.success) ## p.success is prob of registering a single count.

  ## filename = paste("DynNB-synth-", corr.type, "-", P, "-mu-", nb.mean, ".RData", sep="")  
  ## if (FALSE) {
  ##   save(d.true, nb.mean, marg.V, phi.true, W.true, beta, tX, X, log.mean, psi, p.success, y, freq, xgrid,
  ##        file=filename, compress=TRUE)
  ## }

  #-----------------------------------------------------------------------------

  X.dyn = X;

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

  if (est.ar=="with.ar") {
    phi.true = NULL
    W.true = NULL
  }   
  
  out = list()
  
  ## source("Benchmark-DynNB.R")
  for (i in run.idc) {
    ## source("Benchmark-DynNB.R")
    nm = methods[i];
    out[[nm]] <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=X.stc, 
                                  samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                  method=nm, var.names="beta", dset.name="Synth1",
                                  m.0=m0, C.0=C0,
                                  W.a0=W.a0, W.b0=W.b0,
                                  mu.true=mu.true, phi.true=phi.true, W.true=W.true, d.true=NULL);
  }
  
  synth.table = setup.table.dyn(out, "beta")
  
}
