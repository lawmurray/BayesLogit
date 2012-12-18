library("BayesLogit")
library("coda")
library("reglogit")
library("binomlogit")
library("bayesm")

blogit = BayesLogit::logit;

source("LogitPG.R")
source("LogitFS-2010.R")
source("LogitOD.R")
source("Metropolis.R")
source("Benchmark-Utilities.R")

################################################################################
                                  ## SETUP ##
################################################################################

run <- list("synth1"=FALSE,
            "german"=FALSE,
            "ger.num"=FALSE,
            "diabetes"=FALSE,
            "australia"=FALSE,
            "heart"=FALSE,
            "nodal"=FALSE)

write.dir = ""
load.old = FALSE

write.it = FALSE
plot.it  = FALSE
## print.it = TRUE

samp = 10000
burn  = 2000
verbose = 2000
ntrials = 1

logit.meth <- c("PG", "FS", "IndMH", "RAM", "OD",
                "dRUMIndMH", "dRUMHAM", "dRUMAuxMix", "IndivdRUMIndMH", "GP")

run.meth <- logit.meth

start.run = proc.time()

################################################################################
                           ## Dyn Logit Benchmark ##
################################################################################

benchmark.logit <- function(y, X, 
                            samp=1000, burn=100, ntrials=1, verbose=100,
                            method = c("PG", "FS", "GP", "OD", "IndMH",
                              "dRUMIndMH", "dRUMHAM", "dRUMAuxMix", "IndivdRUMIndMH"),
                            m.0=NULL, C.0=NULL,
                            dset.name="", df=Inf)   
{
  
  ## Initialize
  ## sstat.beta.list = list()
  var.names = c("beta");
  ## if (method[1]=="IndMH" || method[1]=="RAM") var.names = c(var.names, "alpha");
  cat("Will benchmark", method[1], "using", dset.name, "dataset for variable(s)", var.names, "\n");

  sstat = list(); for (nm in var.names) sstat[[nm]] = list();
  esstime  = rep(0, ntrials);
  arate    = rep(0, ntrials);
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
      gb$arate = 1
      
    } else if (method=="FS") { ## binary
      
      gb <- logit.mix.gibbs(y, X, samp=samp, burn=burn, b.0=m.0, P.0=P.0, verbose=verbose)
      gb$arate = 1
      
    } else if (method=="GP") { ## binomial
      
      gb <- reglogit(T=sim, y=y, X=X, N=n, flatten=TRUE, sigma=sqrt(diag(C.0)), nu=1,
                     kappa = 1, icept=FALSE, normalize=FALSE, zzero=TRUE,
                     powerprior=FALSE, kmax=442, bstart=NULL, lt=NULL,
                     nup=NULL, method="MH", verb=verbose, burn.point=burn);
      gb$beta = gb$beta[(burn+1):sim,]
      gb$arate = 1
      
    } else if (method=="IndMH") { ## multinomial, fraction
      
      ## This is essentially independent-metropolis in the binary case
      m.0 = array(m.0, dim=c(P, 1));
      P.0 = array(P.0, dim=c(P, P, 1))
      gb <- mlogit.MH.R(y, X, n, m.0, P.0, beta.0=rep(0, P), samp=samp, burn=burn,
                        method="Ind", tune=1.0, df=df, verbose=verbose)
      gb$arate = gb$acceptr
      
    } else if (method=="OD") { ## binary
      
      gb <- logit.OD.R(y, X, m0=m.0, P0=P.0, samp=samp, burn=burn, verbose=verbose)
      gb$arate = 1
      
    } else if (method=="RAM") { ## multinomial
      
      alts = 2
      X.new = createX(alts, na=0, nd=P, Xd=X, Xa=NULL, base=1, INT=FALSE)
      gb <- rmnlIndepMetrop(Data=list(p=alts, y=(y+1), X=X.new),
                            Prior=list(P.0, m.0),
                            Mcmc=list(R=samp, keep=1));
      gb$beta = gb$betadraw;
      gb$alpha = gb$alphadraw;
      gb$ess.time = gb$total.time;
      gb$arate = gb$acceptr
      
    } else if (method=="HH") { ## binary

    } else if (method=="dRUMIndMH") { ## binomial
      
      y = as.numeric(y);
      gb <- dRUMIndMH(y, n, X, sim=sim, burn=burn, b0=m.0, B0=C.0, verbose=verbose);
      gb$beta = t(gb$beta)[(burn+1):sim,];
      gb$ess.time = rep(gb$duration_wBI, 3)
      gb$total.time= gb$duration
      gb$arate = gb$rate / 100;
      
    } else if (method=="dRUMHAM") { ## binomial
      
      y = as.numeric(y)
      gb <- dRUMHAM(y, n, X, sim=sim, burn=burn, b0=m.0, B0=C.0, verbose=verbose);
      gb$beta = t(gb$beta)[(burn+1):sim,];
      gb$ess.time = rep(gb$duration_wBI, 3)
      gb$total.time= gb$duration
      gb$arate = 1
      
    } else if (method=="dRUMAuxMix") { ## binomial - precompute all mixtures
      
      y = as.numeric(y)
      gb <- dRUMAuxMix(y, n, X, sim=sim, burn=burn, b0=m.0, B0=C.0, verbose=verbose);
      gb$beta = t(gb$beta)[(burn+1):sim,];
      gb$ess.time = rep(gb$duration_wBI, 3)
      gb$total.time= gb$duration
      gb$arate = 1
      
    } else if (method=="IndivdRUMIndMH") { ## binary
      
      y = as.numeric(y)
      gb <- IndivdRUMIndMH(y, X, sim=sim, burn=burn, b0=m.0, B0=C.0, verbose=verbose);
      gb$beta = t(gb$beta)[(burn+1):sim,];
      gb$ess.time = rep(gb$duration_wBI, 3)
      gb$total.time= gb$duration
      gb$arate = gb$rate / 100;
      
    } else {
      print("Unknown method.")
      return(NA);
    }
    
    end.time = proc.time();
    gb$call.time = end.time - start.time;

    ## sstat.beta = sum.stat(gb$beta, gb$ess.time[1], thin=1);
    ## sstat.beta.list[[i]] = sstat.beta;

    ## sstat.temp = list();
    for (nm in var.names) { sstat[[nm]][[i]] = sum.stat(gb[[nm]], gb$ess.time[3], thin=1); }

    esstime[i]  = gb$ess.time[3];
    calltime[i] = gb$call.time[3];
    arate[i]    = gb$arate
  }

  for (nm in var.names)
    sstat[[nm]] = simplify2array(sstat[[nm]]);

  ## To get an idea of how ill-conditioned X is.
  eval = eigen(t(X) %*% X)$values;
  cn.X = eval[1] / eval[P]
  
  info <- list("samp" = samp,
               "burn" = burn,
               "ntrials" = ntrials,
               "method" = method,
               "dataset" = dset.name,
               "N" = nrow(X),
               "P" = ncol(X),
               "m.0" = m.0,
               "P.0" = P.0,
               "cn.X"=cn.X,
               "ave.arate"=mean(arate))

  out <- list("gb"=gb, "sstat"=sstat, "ess.time"=esstime, "call.time"=calltime, "info"=info, "arate"=arate);

  out

} ## benchmark.logit

################################################################################

################################################################################

################################################################################
                              ## GEN SYNTHETIC ##
################################################################################

if (FALSE) {

  N = 270;
  P = 14;

  ##------------------------------------------------------------------------------
  ## Correlated predictors
  rho = 0.5
  Sig = matrix(rho, nrow=P, ncol=P); diag(Sig) = 1.0;
  U   = chol(Sig);
  X   = matrix(rnorm(N*P), nrow=N, ncol=P) %*% U;

  ##------------------------------------------------------------------------------
  ## Sparse predictors
  X   = matrix(rnorm(N*P, sd=1), nrow=N, ncol=P);
  vX  = as.numeric(X);
  low  = vX < quantile(vX, 0.5)
  high = vX > quantile(vX, 0.5);
  X[low]  = 0;
  X[!low] = 1;
  
  beta = rnorm(P, mean=0, sd=1);
  
  ## beta = c(1.0, 0.4);
  ## X = matrix(rnorm(N*P), nrow=N, ncol=P);
  
  psi = X %*% beta;
  p = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);

  glm.1 = glm(y ~ X + 0, family=binomial(link=logit))

  save(y, X, N, P, rho, Sig, file="DataSets/synth02-logit.RData")

}

################################################################################
                                 ## TESTING ##
################################################################################

if (FALSE) {

  samp = 10000;
  verbose = 1000;
  burn = 1000;
  ntrials = 1;

  out.pg <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                             method="PG", m.0=NULL, C.0=NULL, dset.name="")

  out.mh <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                            method="IndMH", m.0=NULL, C.0=NULL, dset.name="", df=Inf)

  mh = mlogit.MH.R(y, X, n, m.0, P.0, beta.0=beta.pm, method="Ind", tune=1.0, df=Inf)

  out = list();
  for (i in c(1,5,10)) {
    nm = logit.meth[i]
    out[[nm]] <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                           method=nm, m.0=NULL, C.0=NULL, dset.name="", df=Inf)
  }
  glm.1 = glm(y ~ X + 0, family=binomial(link=logit))

  out.tbl = setup.table(out)
  
  lapply(out, function(x) x$beta)

  lapply(out, function(x) {
    m = min(3, ncol(x$beta));
    if(is.null(x)) return(NULL)
    ave.esr = apply(as.matrix(x$beta[,4,]), 1, mean);
    ss = order(ave.esr)[1:m]
    cbind(x$beta[ss,,], ss)
  });
}

################################################################################
                           ## BENCHMARK DATA SETS ##
################################################################################

##------------------------------------------------------------------------------
## GERMANY
if (run$german) {
  
  load("DataSets/germany.RData");
  dset.name="German"
  file.name="bench-ger.RData"
  
  y = y.ger
  X = X.ger

  glm.ger = glm(y ~ X + 0, family=binomial(link=logit))

  bench.ger = list()
  if (load.old) load(file.name);
  
  out = list();
  for (nm in run.meth) {
    bench.ger[[nm]] <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                                       method=nm, m.0=NULL, C.0=NULL, dset.name=dset.name)
  }

  if (write.it) save(y.ger, X.ger, glm.ger, bench.ger, file=file.name);

}

##------------------------------------------------------------------------------
## GERMAN - NUMERIC
if (run$ger.num) {
  
  load("DataSets/german-numeric.RData");
  dset.name="German.Numeric"
  file.name="bench-ger.num.RData"
  
  y = y.ger.num
  X = X.ger.num

  glm.ger.num = glm(y ~ X + 0, family=binomial(link=logit))

  bench.ger.num = list()
  if (load.old) load(file.name);
  
  for (nm in run.meth) {
    bench.ger.num[[nm]] <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                                           method=nm, m.0=NULL, C.0=NULL, dset.name=dset.name)
  }

  if (write.it) save(y.ger.num, X.ger.num, glm.ger.num, bench.ger.num, file=file.name);

}

##------------------------------------------------------------------------------
## DIABETES
if (run$diabetes) {
  load("DataSets/diabetes.RData");
  dset.name="Diabetes"
  file.name="bench-diabetes.RData"
  
  y = y.diabetes
  X = X.diabetes

  glm.diabetes = glm(y ~ X + 0, family=binomial(link=logit))

  bench.diabetes = list()
  if (load.old) load(file.name);
  
  for (nm in run.meth) {
    bench.diabetes[[nm]] <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                                            method=nm, m.0=NULL, C.0=NULL, dset.name=dset.name)
  }

  if (write.it) save(y.diabetes, X.diabetes, glm.diabetes, bench.diabetes, file=file.name);
  
}

##------------------------------------------------------------------------------
## AUSTRALIA
if (run$australia) {
  
  load("DataSets/australia.RData");
  dset.name="Australia"
  file.name="bench-aus.RData"

  y = y.aus
  X = X.aus

  glm.aus = glm(y ~ X + 0, family=binomial(link=logit))

  bench.aus = list()
  if (load.old) load(file.name);
  
  out = list();
  for (nm in run.meth) {
    bench.aus[[nm]] <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                           method=nm, m.0=NULL, C.0=NULL, dset.name=dset.name)
  }

  if (write.it) save(y.aus, X.aus, glm.aus, bench.aus, file=file.name);
}

##------------------------------------------------------------------------------
## HEART
if (run$heart) {
  
  load("DataSets/heart.RData");
  dset.name = "Heart"
  file.name = "bench-heart.RData"

  y = y.heart
  X = X.heart

  glm.heart = glm(y ~ X + 0, family=binomial(link=logit))

  bench.heart = list()
  if (load.old) load(file.name);
  
  for (nm in run.meth) {
    bench.heart[[nm]] <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                                 method=nm, m.0=NULL, C.0=NULL, dset.name=dset.name)
  }

  if (write.it) save(y.heart, X.heart, glm.heart, bench.heart, file=file.name);
  
}

##------------------------------------------------------------------------------
## NODAL
if (run$nodal) {

  data("nodal")
  dset.name = "Nodal"
  file.name = "bench-nodal.RData"

  y = nodal[, 2]
  X = as.matrix(nodal[,-2])

  y.nodal = y
  X.nodal = X

  glm.nodal = glm(y ~ X + 0, family=binomial(link=logit))

  bench.nodal = list()
  if (load.old) load(file.name);
  
  for (nm in run.meth) {
    bench.nodal[[nm]] <- benchmark.logit(y, X, samp=samp, burn=burn, ntrials=ntrials, verbose=2000,
                                         method=nm, m.0=NULL, C.0=NULL, dset.name=dset.name)
  }

  if (write.it) save(y.nodal, X.nodal, glm.nodal, bench.nodal, file=file.name);
  
}

################################################################################

end.run = proc.time()
run.time = proc.time()

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
