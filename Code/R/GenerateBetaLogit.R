## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


id = c("australia", "germany", "heart", "diabetes");
model.name = "logit"

## Load Data
load("~/RPackage/BayesLogit/Code/R/DataSets/australia.RData")
load("~/RPackage/BayesLogit/Code/R/DataSets/germany.RData")
load("~/RPackage/BayesLogit/Code/R/DataSets/heart.RData")
load("~/RPackage/BayesLogit/Code/R/DataSets/diabetes.RData")

all.data = list(
  aus = list("y"=y.aus, "X"=X.aus),
  ger = list("y"=y.ger, "X"=X.ger),
  heart = list("y"=y.heart, "X"=X.heart),
  diabetes = list("y"=y.diabetes, "X"=X.diabetes)
  )

## Set samples
samp = 10000
burn = 2000
trials = 1

## Set method
methods = c("PG", "FS", "HH", "Pro", "DFAE");
method = "HH"

cat("Method is", method, ".\n");

## Load scipts, etc.
if (method=="PG") {
  library("BayesLogit", lib.loc="~/RPackage/BayesLogit/Code/BLPackage/Test/");
  source("~/RPackage/BayesLogit/Code/R/LogitPG.R");
}
if (method=="FS") source("~/RPackage/BayesLogit/Code/R/LogitFS-2010.R")
if (method=="HH") source("~/RPackage/BayesLogit/Code/R/LogitKS-reduced.R")
if (method=="Pro") {
  library("BayesBridge", lib.loc="~/RPackage/BayesBridge/Code/BBPackage/Test/");
  source("~/RPackage/BayesLogit/Code/R/Probit.R")
}
if (method=="DFAE") source("~/RPackage/BayesLogit/Code/R/LogitScott.R")

## Check how many files we have now.
exst = length( dir(pattern=paste(method, ".*.RData")) )
cat("There are", exst, "similar files.\n");

## Do the simulations...
for (j in 1:4) {

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
      ## out = logit.gibbs(y, X, n.prior=0.0, samp=samp, burn=burn);
      out = logit.gibbs.np.R(y, X, P.0=P.0, samp=samp, burn=burn, verbose=1000);
    } else if (method=="FS") {
      out = logit.mix.gibbs(y, X, P.0=P.0, samp=samp, burn=burn, verbose=1000);
    } else if (method=="HH") {
      out = logit.KS.gibbs(y, X, P.0=P.0, samp=samp, burn=burn, verbose=1000);
    } else if (method=="Pro") {
      out = probit(y, X, samp=samp, burn=burn, verbose=1000, m.0=m.0, C.0=C.0);
    } else if (method=="DFAE") {
      out = logit.DFAE.gibbs(y, X, P.0=P.0, samp=samp, burn=burn, verbose=1000);
    } else {
      print("Bad method.");
      glm1 = glm(y ~ X + 0, family=binomial(link="logit"));
      coef(glm1)
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

  ## I thought maybe using the R compiled package was making things faster.  It
  ## turns out that this is not the case.
  ## dyn.load("~/RPackage/BayesLogit/Code/C/Logit.so")
  ## source("~//RPackage/BayesLogit/Code/C/LogitWrapper.R")
