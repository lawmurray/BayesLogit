library("BayesLogit")
library("coda")

source("NBPG-logmean.R")
source("NBFS-logmean.R")
## source("Metropolis.R")
source("Benchmark-Utilities.R")

################################################################################
                                  ## SETUP ##
################################################################################

run <- list("synth1"=TRUE)

write.dir = ""

write.it = FALSE
plot.it  = FALSE
print.it = TRUE

samp = 10000
burn  = 1000
ntrials = 1

logit.meth <- c("PG", "FS", "IndMH");

################################################################################
                             ## Dyn NB Benchmark ##
################################################################################

benchmark.NB <- function(y, X, 
                            samp=1000, burn=100, ntrials=1, verbose=100,
                            method = c("PG", "FS", "IndMH"),
                            m.0=NULL, C.0=NULL,
                            dset.name="", df=Inf)   
{
  
  ## Initialize
  sstat.beta.list = list();
  esstime  = rep(0, ntrials);
  calltime = rep(0, ntrials);

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

    } else if (method=="FS") { ## NB

      gb <- NB.FS.gibbs(y, X, b.0=m.0, P.0=P.0, samp=samp, burn=burn, verbose=verbose)
      
    } else if (method=="IndMH") { ## multinomial, fraction
      m.0 = array(m.0, dim=c(P, 1));
      P.0 = array(P.0, dim=c(P, P, 1))
    } else {
      print("Unknown method.")
      return(NA);
    }
    
    end.time = proc.time();
    gb$call.time = end.time - start.time;
    
    sstat.beta = sum.stat(gb$beta, gb$ess.time[1], thin=1);

    sstat.beta.list[[i]] = sstat.beta;
    esstime[i]  = gb$ess.time[1];
    calltime[i] = gb$call.time[1];
  }

  sstat.beta = simplify2array(sstat.beta.list);

  info <- list("samp" = samp,
               "burn" = burn,
               "ntrials" = ntrials,
               "method" = method,
               "dataset" = dset.name,
               "N" = nrow(X),
               "P" = ncol(X),
               "m.0" = m.0,
               "P.0" = P.0)

  out <- list("gb"=gb, "beta"=sstat.beta, "esstime"=esstime, "calltime"=calltime, "info"=info);

  out

} ## benchmark.logit

################################################################################
                              ## GEN SYNTHETIC ##
################################################################################

if (FALSE) {

  N = 100;
  P = 2;

  ##------------------------------------------------------------------------------
  ## Highly correlated predictors
  rho = 0.2
  Sig = matrix(rho, nrow=P, ncol=P); diag(Sig) = 1.0;
  U   = chol(Sig);
  X   = matrix(rnorm(N*P), nrow=N, ncol=P) %*% U;

  ##------------------------------------------------------------------------------
  ## Sparse predictors
  X   = matrix(rnorm(N*P, sd=0.0001), nrow=N, ncol=P);
  diag(X) = c(sign(rnorm(P)));
  
  
  beta = rnorm(P, mean=0, sd=2);
  
  ## beta = c(1.0, 0.4);
  ## X = matrix(rnorm(N*P), nrow=N, ncol=P);

  ## log.mean
  mu = exp(X %*% beta);
  d  = 4;
  lambda = (mu / d) * rgamma(N, d, 1);
  y = rpois(N, lambda);

}

################################################################################
                                 ## TESTING ##
################################################################################

if (FALSE) {

  samp = 10000;
  verbose = 1000;
  burn = 1000;
  ntrials = 1;

  out.pg <- benchmark.NB(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                         method="PG", m.0=NULL, C.0=NULL, dset.name="")

  out.fs <- benchmark.NB(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                         method="FS", m.0=NULL, C.0=NULL, dset.name="", df=6)

  out = list();
  for (i in 1:2) {
    out[[i]] <- benchmark.NB(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                           method=logit.meth[i], m.0=NULL, C.0=NULL, dset.name="")
  }

  lapply(out, function(x) x$beta)
  
}
