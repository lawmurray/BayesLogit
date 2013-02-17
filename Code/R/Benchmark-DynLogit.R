## For benchmarking dynamic binomial logistic regression.

## library("BayesLogit")
library("coda")

source("Benchmark-Utilities.R");

source("DynLogitPG.R")
source("DynLogitdRUM.R")
source("CUBS-MH.R")

################################################################################
                                  ## SETUP ##
################################################################################

run <- list("tokyo"=FALSE)

write.dir = ""

write.it = FALSE
plot.it  = FALSE
print.it = FALSE

methods = c("PG", "dRUM", "CUBS", "FS2010");

samp = 1000
burn  = 100
ntrials = 1

################################################################################
                           ## Dyn Logit Benchmark ##
################################################################################

benchmark.dyn.logit <- function(y, X.dyn, n, X.stc=NULL, 
                                samp=1000, burn=100, ntrials=1, verbose=100,
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
      gb <- cubs.mh(y=y, X.dyn=X.dyn, n=rep(n,T), m0=m.0, C0=C.0,
                    samp=samp, burn=burn, verbose=verbose,
                    mu.m0=mu.m0, mu.P0=mu.P0,
                    phi.m0=phi.m0, phi.P0=phi.P0,
                    W.a0=W.a0, W.b0=W.b0, X.stc=X.stc,
                    mu.true = mu.true, phi.true=phi.true, W.true=W.true,
                    obs="binom")
      gb$a.rate = gb$a.rate[samp]
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

    for (nm in var.names) { sstat[[nm]][[i]] = sum.stat.dyn(gb[[nm]], gb$ess.time[3], thin=1); }

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
  
  ## source("Benchmark-DynLogit.R")
  for (i in 1:3) {
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

##------------------------------------------------------------------------------
