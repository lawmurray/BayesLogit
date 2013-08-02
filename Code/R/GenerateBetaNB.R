## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


id = c("synth1");
model.name="NB"

## Load Data
load("~/RPackage/BayesLogit/Code/R/DataSets/NB-synth-1.RData")

all.data = list(
  synth1 = list("y"=NB.synth.1$y, "X"=NB.synth.1$X)
  )
J = length(all.data)

## Set samples
samp = 10000
burn = 2000
trials = 1

## Set method
methods = c("FS", "PG");
method = "FS"

cat("Method is", method, ".\n");

## Load scripts, etc.
if (method=="PG") {
  library("BayesLogit", lib.loc="~/RPackage/BayesLogit/Code/BLPackage/Test/");
  source("~/RPackage/BayesLogit/Code/R/NBPG-logmean.R");
}
if (method=="FS") {
  if(!is.loaded("FSF_nmix.so")) dyn.load("~/RPackage/BayesLogit/Code/R/FSF_nmix.so")
  source("~/RPackage/BayesLogit/Code/R/NBFS-logmean.R")
}

## Check how many files we have now.
exst = length( dir(pattern=paste(method, ".*.RData")) )
cat("There are", exst, "similar files.\n");

## Run simluations...
for (j in 1:J) {

  cat("Working on", id[j], "dataset.\n");
  
  y = all.data[[j]]$y;
  X = all.data[[j]]$X;

  P = ncol(X);
  P.0 = diag(0.01, P);
  C.0 = diag(100, P);
  m.0 = rep(0, P);

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
    
    if (method=="PG") {
      out = NB.PG.gibbs(y, X, P.0=P.0, samp=samp, burn=burn, verbose=1000);
    } else if (method=="FS") {
      out = NB.FS.gibbs(y, X, P.0=P.0, samp=samp, burn=burn, verbose=1000);
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
