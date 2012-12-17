library("BayesLogit")
source("LogitPG.R")
source("Metropolis.R")
source("Logit-MixedModel.R")
source("Benchmark-Utilities.R")

################################################################################

################################################################################

samp = 10000
burn = 1000
verbose = 1000
df = Inf

mm.method = c("PGMM", "IndMHMM", "PG", "IndMH");

var.names <- list("IndMH"=c("ab"),
                  "PG"="ab",
                  "PGMM"=c("ab", "phi"),
                  "IndMHMM"=c("ab", "phi"))

################################################################################

################################################################################

benchmark.blogit.mm <- function(y, X.re, X.fe, n, shape=1, rate=1, m.0=NULL, C.0=NULL,
                                samp=1000, burn=100, ntrials=1, verbose=100,
                                method = c("PGMM", "IndMHMM"),
                                dset.name="", df=Inf, var.names="beta")   
{
  source("Logit-MixedModel.R")
  
  ## Initialize
  cat("Will benchmark", method[1], "using", dset.name, "dataset for variable(s)", var.names, "\n");

  sstat = list(); for (nm in var.names) sstat[[nm]] = list();
  esstime  = rep(0, ntrials);
  arate    = rep(0, ntrials);
  calltime = rep(0, ntrials);

  y = as.matrix(y);
  X.re = as.matrix(X.re);
  X.fe = as.matrix(X.fe);
  X = cbind(X.re, X.fe);
  n = rep(1, length(y));
  N = nrow(X);
  P.a = ncol(X.re)
  P.b = ncol(X.fe)
  P.ab = ncol(X);
  method = method[1];

  ## Use mean precision in a fixed effects model.
  prec.mean = shape / rate;
  m.0.logit = matrix(0, nrow=P.ab, ncol=1)
  P.0.logit = diag(prec.mean, P.ab)
  
  sim = samp + burn;

  gb = list();

  ## Default prior
  if (is.null(m.0)) m.0 = rep(0, P.b);
  if (is.null(C.0)) C.0 = diag(100, P.b);
  P.0 = solve(C.0);
  
  for(i in 1:ntrials) {

    start.time = proc.time();
    
    if (method=="PGMM") { ## Polya Gamma Mixed Model

      gb <- logit.PG.mm(y, X.re, X.fe, n, shape, rate, m.0, P.0, samp=samp, burn=burn, verbose=verbose);
      gb$arate = 1
      
    } else if (method=="IndMHMM") { ## Independence Metrop. Mixed Model.

      gb <- ind.metropolis.blogit(y, X.re, X.fe, n, shape, rate, m.0, P.0, samp=samp, burn=burn, verbose=verbose, df=df)
      gb$arate = gb$acceptr
      
    } else if (method=="PG") { ## binomial, fraction
      
      gb <- logit.R(y, X, n, m0=m.0.logit, P0=P.0.logit, samp=samp, burn=burn, verbose=verbose)
      gb$arate = 1
      gb$ab = gb$beta
      
    } else if (method=="IndMH") { ## multinomial, fraction
      
      ## This is essentially independent-metropolis in the binary case
      m.0.mlogit = array(m.0.logit, dim=c(P.ab, 1));
      P.0.mlogit = array(P.0.logit, dim=c(P.ab, P.ab, 1))
      gb <- mlogit.MH.R(y, X, n, m.0.mlogit, P.0.mlogit, beta.0=rep(0, P.ab), samp=samp, burn=burn,
                        method="Ind", tune=1.0, df=df, verbose=verbose)
      gb$arate = gb$acceptr
      gb$ab = gb$beta
      
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

  eval = eigen(t(X) %*% X)$values;
  cn.X = eval[1] / eval[P.ab]
  
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

if (FALSE) {

  load("DataSets/mm-example.RData")
  dset.name = "mm-example"
  
  P.a = ncol(X.re);
  P.b = ncol(X.fe);

  shape = 1
  rate  = 1
  
  prec.b = rep(0.01, P.b)
  
  m.0   = rep(0, P.b);
  P.0   = diag(prec.b, P.b);

  samp = 10000
  burn = 1000
  verbose = 1000
  df = Inf
  ntrials = 1
  
  out = list();
  
  for (i in 1:4) {
    nm = mm.method[i]
    out[[nm]] <- benchmark.blogit.mm(y, X.re, X.fe, n, shape, rate, m.0=m.0, C.0=solve(P.0),
                                     samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                     method=nm,
                                     dset.name=dset.name, df=df, var.names=var.names[[nm]])
  }

  setup.table(out, "ab")
  
}
