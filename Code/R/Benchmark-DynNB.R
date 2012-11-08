library("BayesLogit")
library("coda")

source("Benchmark-Utilities.R")

source("DynNBPG.R")
source("DynNBFS-2009.R")

################################################################################
                                  ## SETUP ##
################################################################################

run <- list("flu"   =FALSE,
            "synth1"=FALSE)

write.dir = ""

write.it = FALSE
plot.it  = FALSE
print.it = FALSE

samp = 100
burn  = 100
ntrials = 1

################################################################################
                             ## Dyn NB Benchmark ##
################################################################################

benchmark.dyn.NB <- function(y, X.dyn, X.stc=NULL,
                             samp=1000, burn=100, ntrials=1, verbose=100,
                             method="PG",
                             m.0=NULL, C.0=NULL,
                             mu.m0=NULL, mu.V0=NULL,
                             phi.m0=NULL, phi.V0=NULL,
                             W.a0=NULL, W.b0=NULL,
                             d.true=NULL,
                             beta.true=NULL, iota.true=NULL,
                             mu.true=NULL, phi.true=NULL, W.true=NULL)
{

  sstat.beta.list = list();
  sstat.iota.list = list();

  info.beta = matrix(nrow=ntrials, ncol=5);
  colnames(info.beta) = c("ave.ESS", "sd.ESS", "ave.ESS.sec", "sd.ESS.sec", "ess.time");
  info.iota = info.beta;

  for(i in 1:ntrials) {

    if (method=="PG") {
      gb <- dyn.NB.PG(y=y, X.dyn=X.dyn, X.stc=X.stc,
                      samp=samp, burn=burn, verbose=verbose,
                      m.0=m.0, C.0=C.0,
                      mu.m0=mu.m0, mu.V0=mu.V0,
                      phi.m0=phi.m0, phi.V0=phi.V0,
                      W.a0=W.a0, W.b0=W.b0,
                      d.true=d.true, w.true=NULL,
                      beta.true=beta.true, iota.true=iota.true,
                      mu.true=mu.true, phi.true=phi.true, W.true=W.true)
    } else if (method=="FS") {
      gb <- dyn.NB.FS(y=y, X.dyn=X.dyn, X.stc=X.stc,
                      samp=samp, burn=burn, verbose=verbose,
                      m.0=m.0, C.0=C.0,
                      mu.m0=mu.m0, mu.V0=mu.V0,
                      phi.m0=phi.m0, phi.V0=phi.V0,
                      W.a0=W.a0, W.b0=W.b0,
                      d.true=d.true, lambda.true=NULL, r.true=NULL,
                      beta.true=beta.true, iota.true=iota.true,
                      mu.true=mu.true, phi.true=phi.true, W.true=W.true)
    }
    else {
      print("Unknown method.")
      return(NA);
    }

    sstat.beta = sum.stat.dyn(gb$beta, gb$ess.time[1], thin=1);
    sstat.iota = sum.stat    (gb$iota, gb$ess.time[1], thin=1);

    info.beta[i,"ave.ESS"]     = mean(sstat.beta[,,"ESS"]);
    info.beta[i,"sd.ESS" ]     = sd(  sstat.beta[,,"ESS"]);
    info.beta[i,"ave.ESS.sec"] = mean(sstat.beta[,,"ESS.sec"]);
    info.beta[i,"sd.ESS.sec"]  = sd(  sstat.beta[,,"ESS.sec"]);
    info.beta[i,"ess.time"]    = gb$ess.time[1]

    info.iota[i,"ave.ESS"]     = mean(sstat.iota[,"ESS"]);
    info.iota[i,"sd.ESS" ]     = sd(  sstat.iota[,"ESS"]);
    info.iota[i,"ave.ESS.sec"] = mean(sstat.iota[,"ESS.sec"]);
    info.iota[i,"sd.ESS.sec"]  = sd(  sstat.iota[,"ESS.sec"]);
    info.iota[i,"ess.time"]    = gb$ess.time[1]

    sstat.beta.list[[i]] = sstat.beta;
    sstat.iota.list[[i]] = sstat.iota;

  }

  sstat.beta = simplify2array(sstat.beta.list);
  sstat.iota = simplify2array(sstat.iota.list);

  out <- list("gb"=gb, "sstat.beta"=sstat.beta, "info.beta"=info.beta,
              "sstat.iota"=sstat.iota, "info.iota"=info.iota);

  out

} ## benchmark.dyn.NB

################################################################################
                           ## BENCHMARK DATA SETS ##
################################################################################

##------------------------------------------------------------------------------
                                 ## FLU DATA ##
##------------------------------------------------------------------------------

if (run$flu) {

}

##------------------------------------------------------------------------------
                               ## SYNTHETIC 1 ##
##------------------------------------------------------------------------------

if (run$synth1) {

  load("Benchmark-DataSets/DynNB-synth-1.RData")

  y = dyn.nb.1$y;
  X.dyn = dyn.nb.1$X;

  ## Prior
  b.m0 = 3.0;
  b.C0 = 3.0;
  W    = 0.1;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  ## source("Benchmark-DynNB.R")
  pg <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=NULL,
                         samp=samp, burn=burn, ntrials=ntrials, verbose=100,
                         method="PG",
                         m.0=b.m0, C.0=b.C0,
                         W.a0=W.a0, W.b0=W.b0,
                         mu.true=0.0, phi.true=1.0)

  ## source("Benchmark-DynNB.R")
  fs <- benchmark.dyn.NB(y, X.dyn=X.dyn, X.stc=NULL,
                         samp=samp, burn=burn, ntrials=ntrials, verbose=100,
                         method="FS",
                         m.0=b.m0, C.0=b.C0,
                         W.a0=W.a0, W.b0=W.b0,
                         mu.true=0.0, phi.true=1.0)

  if (print.it) { print("PG beta: "); print(pg$info.beta); print("FS beta:"); print(fs$info.beta); }
  if (plot.it)  { plot.bench(pg, fs); plot.check.NB(y, X.dyn, bmark1=pg, bmark2=fs); }
  if (write.it) save(pg, fs, file=file.path(write.dir, "bmark-synth1.RData"))
  
}


##------------------------------------------------------------------------------
