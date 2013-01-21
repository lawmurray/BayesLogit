library("BayesLogit")
library("coda")
library("bayesm")

source("NBPG-logmean.R")
source("NBFS-logmean.R")
## source("Metropolis.R")
source("Benchmark-Utilities.R")

################################################################################

################################################################################

################################################################################
                                  ## SETUP ##
################################################################################

run <- list("small"=FALSE,
            "med"=FALSE)

write.dir = ""

write.it = FALSE
plot.it  = FALSE
print.it = TRUE

samp = 10000
burn  = 2000
ntrials = 10
verbose = 1000

logit.meth <- c("PG", "FS", "RAM");

################################################################################
                             ## Dyn NB Benchmark ##
################################################################################

benchmark.NB <- function(y, X, 
                         samp=1000, burn=100, ntrials=1, verbose=100,
                         method = c("PG", "FS", "IndMH"),
                         m.0=NULL, C.0=NULL,
                         dset.name="")   
{
  
  ## Initialize
  var.names="beta"
  sstat = list(); for (nm in var.names) sstat[[nm]] = list();
  ess.time  = rep(0, ntrials);
  call.time = rep(0, ntrials);
  arate     = rep(0, ntrials);
  
  y = as.matrix(y);
  X = as.matrix(X);
  n = rep(1, length(y));
  N = nrow(X);
  P = ncol(X);
  method = method[1];

  sim = samp + burn;

  gb = list();

  ## Default prior
  if (is.null(m.0)) m.0 = rep(0, P);
  if (is.null(C.0)) C.0 = diag(100, P);
  P.0 = solve(C.0);
  
  for(i in 1:ntrials) {

    start.time = proc.time();
    
    if (method=="PG") { ## NB

      gb <- NB.PG.gibbs(y, X, b.0=m.0, P.0=P.0, samp=samp, burn=burn, verbose=verbose)
      gb$arate = 1

    } else if (method=="FS") { ## NB

      gb <- NB.FS.gibbs(y, X, b.0=m.0, P.0=P.0, samp=samp, burn=burn, verbose=verbose)
      gb$arate = 1
      
    } else if (method=="RAM") {

      gb <- rnegbinRw(Data=list(y=y,X=X), Prior=list(betabar=m.0, A=P.0, a=0.001, b=0.001), Mcmc=list(R=sim, keep=1, burn=burn));
      gb$beta = gb$betadraw[(burn+1):sim,]
      gb$total.time = proc.time() - start.time;
      gb$arate = gb$acceptrbeta
      
    } else {
      print("Unknown method.")
      return(NA);
    }
    
    end.time = proc.time();
    gb$call.time = end.time - start.time;

    sstat.temp = list();
    for (nm in var.names) { sstat[[nm]][[i]] = sum.stat(gb[[nm]], gb$ess.time[3], thin=1); }

    ess.time[i]  = gb$ess.time[3];
    call.time[i] = gb$call.time[3];
    arate[i]     = gb$arate
  }

  for (nm in var.names)
    sstat[[nm]] = simplify2array(sstat[[nm]]);
  
  info <- list("samp" = samp,
               "burn" = burn,
               "ntrials" = ntrials,
               "method" = method,
               "dataset" = dset.name,
               "N" = nrow(X),
               "P" = ncol(X),
               "m.0" = m.0,
               "P.0" = P.0,
               "ave.arate"=mean(arate))

  out <- list("gb"=gb, "sstat"=sstat, "ess.time"=ess.time, "call.time"=call.time, "info"=info);

  out

} ## benchmark.logit

################################################################################
                              ## GEN SYNTHETIC ##
################################################################################

if (FALSE) {

  N = 400;
  P = 4;

  ##------------------------------------------------------------------------------
  ## Correlated predictors
  rho = 0.2
  Sig = matrix(rho, nrow=P, ncol=P); diag(Sig) = 1.0;
  U   = chol(Sig);
  X   = matrix(rnorm(N*P), nrow=N, ncol=P) %*% U;

  ##------------------------------------------------------------------------------
  ## Sparse predictors
  ## X   = matrix(rnorm(N*P, sd=0.0001), nrow=N, ncol=P);
  ## diag(X) = c(sign(rnorm(P)));

  ## Use an intercpet
  X[,1] = 1.0
  icept = 3.0
  beta = rnorm(P, mean=0, sd=icept / 10);
  beta[1] = icept
  
  ## beta = c(1.0, 0.4);
  ## X = matrix(rnorm(N*P), nrow=N, ncol=P);

  ## log.mean
  mu = exp(X %*% beta);
  d  = 4;
  lambda = (mu / d) * rgamma(N, d, 1);
  y = rpois(N, lambda);

  ## save(N, P, rho, Sig, X, mu, d, lambda, y, beta, file="NBSynth_N400_P4_I3.RData")

}

################################################################################
                                 ## TESTING ##
################################################################################

if (FALSE) {

  samp = 10000;
  verbose = 2000;
  burn = 1000;
  ntrials = 10;

  out.pg <- benchmark.NB(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                         method="PG", m.0=NULL, C.0=NULL, dset.name="")

  out.fs <- benchmark.NB(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                         method="FS", m.0=NULL, C.0=NULL, dset.name="")

  out = list();
  for (i in 1:3) {
    nm = logit.meth[i]
    out[[nm]] <- benchmark.NB(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                           method=nm, m.0=NULL, C.0=NULL, dset.name="")
  }

  setup.table(out)
  
}

################################################################################

################################################################################

##------------------------------------------------------------------------------
## SMALL COUNTS
if (run$small)
{
  load("DataSets/NBSynth_N400_P4_I2.RData")

  out = list();
  for (i in 1:3) {
    nm = logit.meth[i]
    out[[nm]] <- benchmark.NB(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                           method=nm, m.0=NULL, C.0=NULL, dset.name="")
  }

  bench.nb.small.counts = out

  if (write.it) save(bench.nb.small.counts, file="NB-bench-I2.RData")
}

##------------------------------------------------------------------------------
## MEDIUM COUNTS
if (run$med)
{
  load("DataSets/NBSynth_N400_P4_I3.RData")

  out = list();
  for (i in 1:3) {
    nm = logit.meth[i]
    out[[nm]] <- benchmark.NB(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                           method=nm, m.0=NULL, C.0=NULL, dset.name="")
  }

  bench.nb.med.counts = out

  if (write.it) save(bench.nb.small.counts, file="NB-bench-I3.RData")
}

################################################################################

################################################################################

if (FALSE) {

  setup.table(out)
  
}
