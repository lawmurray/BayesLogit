## library("BayesLogit")
library("coda")

source("Benchmark-Utilities.R")

source("DynNBPG.R")
source("DynNBFS-2009.R")
source("DynNBCUBS.R")

################################################################################
                                  ## SETUP ##
################################################################################

run <- list("flu"   =FALSE,
            "synth1"=FALSE,
            "synth2"=FALSE)

methods = c("PG", "FS", "CUBS")

write.dir = ""

write.it = FALSE
plot.it  = FALSE

samp = 1000
burn  = 100
verbose = 100
ntrials = 1

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
                             mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## n is total number of trials per observation - not vectorized, i.e. the same for each i.
  cat("Will benchmark", method[1], "using", dset.name, "dataset for variable(s)", var.names, "\n");

  sstat = list(); for (nm in var.names) sstat[[nm]] = list();
  arate    = rep(0, ntrials);
  ess.time = rep(0, ntrials);
  
  for(i in 1:ntrials) {

    if (method=="PG") {
      gb <- dyn.NB.PG(y=y, X.dyn=X.dyn, X.stc=X.stc,
                      samp=samp, burn=burn, verbose=verbose,
                      m.0=m.0, C.0=C.0,
                      mu.m0=mu.m0, mu.P0=mu.P0,
                      phi.m0=phi.m0, phi.P0=phi.P0,
                      W.a0=W.a0, W.b0=W.b0,
                      d.true=d.true, w.true=NULL,
                      beta.true=beta.true, iota.true=iota.true,
                      mu.true=mu.true, phi.true=phi.true, W.true=W.true)
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
      gb$a.rate=1
    } else if (method=="CUBS") {
      gb <- dyn.NB.CUBS(y=y, X.dyn=X.dyn, m0=m.0, C0=C.0,
                        samp=samp, burn=burn, verbose=verbose,
                        mu.m0=mu.m0, mu.P0=mu.P0,
                        phi.m0=phi.m0, phi.P0=phi.P0,
                        W.a0=W.a0, W.b0=W.b0, X.stc=X.stc,
                        mu.true = mu.true, phi.true=phi.true, W.true=W.true, d.true=d.true)
      gb$a.rate = gb$a.rate[samp]
    } else {
      print("Unknown method.")
      return(NA);
    }

    for (nm in var.names) { sstat[[nm]][[i]] = sum.stat.dyn(gb[[nm]], gb$ess.time[3], thin=1); }

    arate[i]    = gb$a.rate
    ess.time[i] = gb$ess.time[3]

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

  ## Prior
  b.m0 = 3.0;
  b.C0 = 3.0;
  W    = 0.1;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  bench.synth1 = list();
  
  ## source("Benchmark-DynNB.R")
  for (i in 1:3) {
    ## source("Benchmark-DynNB.R")
    nm = methods[i];
    bench.synth1[[nm]] <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=NULL, 
                                     samp=samp, burn=burn, ntrials=1, verbose=verbose,
                                     method=nm, var.names="beta", dset.name="Synth1",
                                     m.0=b.m0, C.0=b.C0,
                                     W.a0=W.a0, W.b0=W.b0,
                                     mu.true=0.0, phi.true=1.0, d.true=NULL);
  }
  
  synth1.table = setup.table.dyn(bench.synth1, "beta")
 
  if (plot.it)  { plot.bench(pg, fs); plot.check.NB(y, X.dyn, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.synth1, synth1.table, file=file.path(write.dir, "bmark-synth1.RData"))
  
}

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
  for (i in 1:3) {
    ## source("Benchmark-DynNB.R")
    nm = methods[i];
    bench.synth2[[nm]] <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=NULL, 
                                     samp=samp, burn=burn, ntrials=1, verbose=verbose,
                                     method=nm, var.names="beta", dset.name="Synth2",
                                     m.0=b.m0, C.0=b.C0,
                                     W.a0=W.a0, W.b0=W.b0,
                                     mu.true=0.0, phi.true=1.0);
  }
  
  synth2.table = setup.table.dyn(bench.synth2, "beta")
 
  if (plot.it)  { plot.bench(pg, fs); plot.check.NB(y, X.dyn, bmark1=pg, bmark2=fs); }
  if (write.it) save(bench.synth2, synth2.table, file=file.path(write.dir, "bmark-synth2.RData"))

}

##------------------------------------------------------------------------------
