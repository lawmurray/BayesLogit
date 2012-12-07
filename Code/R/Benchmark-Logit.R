library("BayesLogit")
library("coda")
library("reglogit")
library("binomlogit")
library("bayesm")

blogit = BayesLogit::logit;

source("LogitPG.R")
source("LogitFS-2010.R")
source("Metropolis.R")
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

logit.meth <- c("PG", "FS", "GP", "IndMH",
                "dRUMIndMH", "dRUMHAM", "dRUMAuxMix", "IndivdRUMIndMH")

################################################################################
                             ## Dyn NB Benchmark ##
################################################################################

benchmark.logit <- function(y, X, 
                            samp=1000, burn=100, ntrials=1, verbose=100,
                            method = c("PG", "FS", "GP", "OD", "IndMH",
                              "dRUMIndMH", "dRUMHAM", "dRUMAuxMix", "IndivdRUMIndMH"),
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
    
    if (method=="PG") { ## binomial, fraction
      gb <- logit.R(y, X, n, m0=m.0, P0=P.0, samp=samp, burn=burn, verbose=verbose)
    } else if (method=="FS") { ## binary
      gb <- logit.mix.gibbs(y, X, samp=samp, burn=burn, b.0=m.0, P.0=P.0, verbose=verbose)
    } else if (method=="GP") { ## binomial
      gb <- reglogit(T=samp, y=y, X=X, N=n, flatten=TRUE, sigma=1, nu=1,
                     kappa = 0, icept=FALSE, normalize=FALSE, zzero=TRUE,
                     powerprior=FALSE, kmax=442, bstart=NULL, lt=NULL,
                     nup=NULL, method="MH", verb=verbose);
      gb$ess.time = NA
    } else if (method=="OD") { ## binary
      return(NA)
      ## This is essentially independent-metropolis in the binary case
    } else if (method=="IndMH") { ## multinomial, fraction
      m.0 = array(m.0, dim=c(P, 1));
      P.0 = array(P.0, dim=c(P, P, 1))
      gb <- mlogit.MH.R(y, X, n, m.0, P.0, beta.0=rep(0, P), method="Ind", tune=1.0, df=df)
    } else if (method=="RAM") { ## multinomial
      return(NA)
      ## bins = 2
      ## X.new = matrix(0, nrow = N * bins, ncol = bins * P);
      ## for (i in 1:bins) X.new[i+bins*(0:(N-1)),1:P+(i-1)*P] = X;      
      ## gb <- rmnlIndepMetrop(Data=list(p=bins, y=(y+1), X=X.new),
      ##                       Prior=list(P.0, m.0),
      ##                       Mcmc=list(R=sim, keep=1));
      ## gb$beta = gb$betadraw[(burn+1):sim,];
    } else if (method=="HH") { ## binary

    } else if (method=="dRUMIndMH") { ## binomial
      y = as.numeric(y);
      gb <- dRUMIndMH(y, n, X, sim=sim, burn=burn, b0=m.0, B0=C.0, verbose=verbose);
      gb$beta = t(gb$beta)[(burn+1):sim,];
      gb$ess.time = gb$duration_wBI
      gb$total.time= gb$duration
    } else if (method=="dRUMHAM") { ## binomial
      y = as.numeric(y)
      gb <- dRUMHAM(y, n, X, sim=sim, burn=burn, b0=m.0, B0=C.0, verbose=verbose);
      gb$beta = t(gb$beta)[(burn+1):sim,];
      gb$ess.time = gb$duration_wBI
      gb$total.time= gb$duration
    } else if (method=="dRUMAuxMix") { ## binomial - precompute all mixtures
      y = as.numeric(y)
      gb <- dRUMAuxMix(y, n, X, sim=sim, burn=burn, b0=m.0, B0=C.0, verbose=verbose);
      gb$beta = t(gb$beta)[(burn+1):sim,];
      gb$ess.time = gb$duration_wBI
      gb$total.time= gb$duration
    } else if (method=="IndivdRUMIndMH") { ## binary
      y = as.numeric(y)
      gb <- IndivdRUMIndMH(y, X, sim=sim, burn=burn, b0=m.0, B0=C.0, verbose=verbose);
      gb$beta = t(gb$beta)[(burn+1):sim,];
      gb$ess.time = gb$duration_wBI
      gb$total.time= gb$duration
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

  N = 40;
  P = 20;

  ##------------------------------------------------------------------------------
  ## Highly correlated predictors
  rho = 0.99
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
  
  psi = X %*% beta;
  p = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);

}

################################################################################
                                 ## TESTING ##
################################################################################

if (FALSE) {

  N = 200;
  
  beta = c(1.0, 0.4);
  X = cbind(1, rnorm(N));

  psi = X %*% beta;
  p = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);

  samp = 10000;
  verbose = 1000;
  burn = 1000;
  ntrials = 1;

  out.pg <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                             method="PG", m.0=NULL, C.0=NULL, dset.name="")

  out.mh <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                         method="IndMH", m.0=NULL, C.0=NULL, dset.name="", df=6)

  out = list();
  for (i in 1:7) {
    out[[i]] <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                           method=logit.meth[i], m.0=NULL, C.0=NULL, dset.name="")
  }

  lapply(out, function(x) x$beta)
  
}

################################################################################
                           ## BENCHMARK DATA SETS ##
################################################################################



################################################################################
                                 ## APPENDIX ##
################################################################################

## All of the Fussl code is in R.  They precompute the mixture of normals.  If
## we can beat them in the binary case, then we should be able to beat them in
## binomial.

## dRUM stands for difference random utility model.

## dRUMHAM uses a normal apprixmation (like in dRUMIndMH) for values of y_i/n_i
## that are not in the tail and to use the discrete mixture of normals when y_i
## / n_i takes on relatively extreme values.

## dRUMAuxMix uses a discrete mixture of normals.

## dRUMIndMH approximates the error term from the dRUM representation using a
## single normal random variable (matched by mean and variance) as a proposal in
## a Ind. MH sampler.

## IndivdRUMIndMH is like the above for binary logistic.

## In both of the MH steps, you still use data augmentation, but you approximate
## the error terms using a normal.  In this way it is still following the
## discrete mixture of normals approach, but the discrete mixture has only a
## single term.  Unlike the discrete mixture approach, they actually adjust
## using MH, whereas before, they claim that the approximation works well enough
## without it.
