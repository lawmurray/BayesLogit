id = c("synth1");
model.name = "DynLogit"

## Load data
load("~/RPackage/BayesLogit/Code/R/DataSets/DynLogit-synth-1.RData")

all.data = list(
  synth1 = list("y"=dyn.logit.1$y, "X"=dyn.logit.1$X, "with.iota"=dyn.logit.1$with.iota)
  )
J = length(all.data)

## Set samples
samp = 100
burn = 20
trials = 1

## Set method
methods = c("FS", "PG");
method = "FS"

cat("Method is", method, ".\n");

## Load scripts, etc.
if (method=="PG") {
  library("BayesLogit", lib.loc="~/RPackage/BayesLogit/Code/BLPackage/Test/");
  source("~/RPackage/BayesLogit/Code/R/DynLogitPG.R");
}
if (method=="FS") {
  if(!is.loaded("~/RPackage/BayesLogit/Code/R/FSF_nmix.so"))
    dyn.load("~/RPackage/BayesLogit/Code/R/FSF_nmix.so")
  source("~/RPackage/BayesLogit/Code/R/DynBinaryLogitFS-2010.R")
}

## Check how many files we have now.
exst = length( dir(pattern=paste(method, ".*.RData")) )
cat("There are", exst, "similar files.\n");

## Run simulations...
for (j in 1:J) {

  cat("Working on", id[j], "dataset.\n");
  
  y = all.data[[j]]$y;
  X = all.data[[j]]$X;
  with.iota = all.data[[j]]$with.iota

  P = ncol(X);
  
  ## Prior 
  m.0    = rep(0, P+with.iota);
  C.0    = diag(9, P+with.iota);
  W.a0   = rep(10, P);
  W.b0   = rep(1.0, P);
  iota   = NULL; if(!with.iota) iota=0;

  ## Known
  mu  = rep(0, P)
  phi = rep(1, P)

  if (!exists("the.beta")) {
    the.beta = list();
    the.time = list();
    ess.time = list();
    the.mem  = list();
  }

  for (i in 1:trials) {
    cat("Working on trial", i, "using", method, "\n");
    ## Give all the same prior on beta...
    ## Just comment out the two types or routines you do not want to run.

    file.name = paste(model.name, method, id[j], "RData", sep=".");
    
    if (method=="PG") { ## POLYA GAMMA ##
      out <- dyn.logit.PG(y, X, samp=samp, burn=burn, verbose=100,
                          m.0=m.0, C.0=C.0,
                          mu.m0=NULL, mu.V0=NULL,
                          phi.m0=NULL, phi.V0=NULL,
                          W.a0=W.a0, W.b0=W.b0,
                          beta.true=NULL, iota.true=iota, w.true=NULL,
                          mu.true=mu, phi.true=phi, W.true=NULL)
    } else if (method=="FS") { ## FRUHWIRTH-SCHNATTER ##
      out <-dyn.logit.FS(y, X,
                         samp=samp, burn=burn, verbose=100,
                         m.0=m.0, C.0=C.0,
                         mu.m0=NULL, mu.V0=NULL,
                         phi.m0=NULL, phi.V0=NULL,
                         W.a0=W.a0, W.b0=W.b0,
                         z.true=NULL, r.true=NULL,
                         beta.true=NULL, iota.true=iota,
                         mu.true=0.0, phi.true=1.0, W.true=NULL)
    } else {
      print("Bad method.");
      ## glm1 = glm(y ~ X + 0, family=binomial(link="logit"));
      ## coef(glm1)
      break;
    }

    the.beta[[i]] = out$beta;
    the.time[[i]] = out$total.time;
    ess.time[[i]] = out$ess.time
    the.mem[[i]]  = object.size(out);

    ## If using RunGenData.bash
    ## file.name=paste(i, file.name, sep=".");
    ## beta = out$beta;
    ## save(beta, total.time, file=file.name, compress=TRUE);

  }

  all.mem = sapply(ls(), function(x){object.size(get(x))})
  ## For KB, MB, GB divide by 2^10, 2^20, 2^30.
  
  ## If not using RunGenData.bash
  ##saveall save(the.beta, the.time, ess.time, the.mem, all.mem, file=file.name, compress=TRUE);
  
}

################################################################################
                                 ## APPENDIX ##
################################################################################
