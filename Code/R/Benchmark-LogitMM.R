library("BayesLogit")
library("lme4")
library("blme")
source("LogitPG.R")
source("Metropolis.R")
source("Logit-MixedModel.R")
source("Benchmark-Utilities.R")

################################################################################
                                  ## SET UP ##
################################################################################

run <- list("ex01"= FALSE,
            "polls" = FALSE,
            "xerop" = FALSE)

samp = 1000
burn = 200
verbose = 100
kappa = Inf
df = Inf
ntrials = 1

mm.method = c("PGMM", "IndMM1", "IndMM2", "PGMM2", "IndMM3", "IndMM4", "PG", "IndMH");
run.methods = 1:7;

write.it = FALSE
load.old = FALSE

var.names <- list("PGMM"=c("ab", "phi"),
                  "IndMM1"=c("ab", "phi"),
                  "IndMM2"=c("ab", "phi"),
                  "IndMM3"=c("ab", "theta"),
                  "IndMM4"=c("ab", "phi"),
                  "PG"="ab",
                  "IndMH"="ab")

################################################################################
                            ## BENCHMARK ROUTINE ##
################################################################################

benchmark.blogit.mm <- function(y, X.re, X.fe, n, shape=1, rate=1, m.0=NULL, C.0=NULL,
                                samp=1000, burn=100, ntrials=1, verbose=100,
                                method = c("PGMM", "IndMM1"),
                                dset.name="", df=Inf, var.names="beta", kappa=0)   
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
      
    } else if (method=="PGMM2") { ## Polya Gamma Mixed Model

      gb <- logit.PG.mm.2(y, X.re, X.fe, n, shape, rate, kappa=kappa, m.0, P.0, samp=samp, burn=burn, verbose=verbose);
      gb$arate = 1
      
    } else if (method=="IndMM1") { ## Independence Metrop. Mixed Model.

      gb <- ind.metropolis.blogit(y, X.re, X.fe, n, shape, rate, m.0, P.0, samp=samp, burn=burn, verbose=verbose, df=df)
      gb$arate = gb$acceptr
      
    } else if (method=="IndMM2") {

      gb <- ind.metropolis.blogit.2(y, X.re, X.fe, n, shape, rate, m.0, P.0, samp=samp, burn=burn, verbose=verbose, df=df)
      gb$arate = gb$acceptr
      
    } else if (method=="IndMM3") {

      gb <- ind.metropolis.blogit.3(y, X.re, X.fe, n, shape, rate, kappa=kappa, m.0, P.0,
                                    samp=samp, burn=burn, verbose=verbose, df=df)
      gb$arate = gb$acceptr
      
    } else if (method=="IndMM4") {

      gb <- ind.metropolis.blogit.4(y, X.re, X.fe, n, shape, rate, kappa=kappa, m.0, P.0,
                                    samp=samp, burn=burn, verbose=verbose, df=df)
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
                                ## BENCHMARKS ##
################################################################################

##------------------------------------------------------------------------------
## SYNTHETIC DATA
if (run$ex01) {

  load("DataSets/mm-example.RData")
  dset.name = "mm-example"
  file.name = "bench-logitmm-ex01.RData"
  
  P.a = ncol(X.re);
  P.b = ncol(X.fe);

  shape = 2
  rate  = 2
  
  prec.b = rep(0.01, P.b)
  
  m.0   = rep(0, P.b);
  P.0   = diag(prec.b, P.b);

  ## samp = 10000
  ## burn = 2000
  ## verbose = 1000
  ## df = Inf
  ## ntrials = 1
  
  bench.ex01 = list()
  if (load.old) load(file.name);
  
  for (i in run.methods) {
    nm = mm.method[i]
    bench.ex01[[nm]] <- benchmark.blogit.mm(y, X.re, X.fe, n, shape, rate, m.0=m.0, C.0=solve(P.0),
                                            samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                            method=nm,
                                            dset.name=dset.name, df=df, var.names=var.names[[nm]], kappa=kappa)
  }
  
  ex01.tbl = setup.table(bench.ex01, "ab")
  
  if (write.it) save(y, X.re, X.fe, n, shape, rate, prec.b, m.0, P.0, df, bench.ex01, ex01.tbl, file=file.name);

  ## .........................................................
  ## Check sensitivity to hyperparameters and model structure.

  kappa = 1.0; maxit=1000
  ## Without mean of random effects.
  check <- blogit.mm.map(y, X.re, X.fe, n, shape=shape, rate=rate, 
                         m.0=m.0, P.0=P.0, abphi.0=NULL, calc.V=TRUE, trace=FALSE, maxit=maxit)
  c(check$m[1:5], check$m[6], check$m[7])
  
  ## With mean of random effects. log gamma prior
  check <- blogit.mm.map.3(y, X.re, X.fe, n, shape=shape, rate=rate, kappa=kappa,
                           m.0=m.0, P.0=P.0, abtm.0=NULL, calc.V=TRUE, trace=FALSE, maxit=maxit)
  c(check$m[1:5], check$m[6], exp(check$m[7]), check$m[8])
  blogit.llh.mm.3(check$m, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0)

  abtm.0 = NULL
  ## abtm.0 = check$m; abtm.0[7] = exp(abtm.0[7])
  ## With mean of random effects.  gamma prior.
  ## shape = 2, rate=2, kappa=100 produces result almost like hglm1 when using Nelder-Mead
  ## This is not the case when using CG.
  check <- blogit.mm.map.4(y, X.re, X.fe, n, shape=shape, rate=rate, kappa=kappa,
                           m.0=m.0, P.0=P.0, abphim.0=NULL, calc.V=TRUE, trace=FALSE, maxit=maxit)
  c(check$m[1:5], check$m[6], check$m[7], check$m[8])
  blogit.llh.mm.4(check$m, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0)

  abtm.0 = NULL
  ## abtm.0 = check$m; abtm.0[7] = log(abtm.0[7])
  check <- blogit.mm.map.5(y, X.re, X.fe, n, shape=shape, rate=rate, kappa=kappa,
                           m.0=m.0, P.0=P.0, abtm.0=abtm.0, calc.V=TRUE, trace=FALSE, maxit=maxit)
  c(check$m[1:5], check$m[6], exp(check$m[7]), check$m[8])
  blogit.llh.mm.5(check$m, y, X.re, X.fe, n, shape, rate, kappa, m.0, P.0)
  
  ## Check results using lme4 package.  I'm still unclear on their priors.
  mm.df = data.frame("y"=y, "FE1"=X.fe[,1], "group"=factor(rowSums(t(t(X.re) * 1:5))))
  hglm1 = glmer(y ~ FE1 + (1 | group), family=binomial(link="logit"), data=mm.df)
  hglm2 = glmer(y ~ FE1 + (1 | group), family=binomial(link="logit"), data=mm.df, REML=FALSE)
  coef(hglm1)
  coef(hglm2)
  r = ranef(hglm1, postVar = TRUE)
  summary(hglm1)
  attr(r$group,"postVar")

  blme1 = bglmer(y ~ FE1 + (1 | group), family=binomial(link="logit"), data=mm.df)
  
}

##------------------------------------------------------------------------------
## POLLING DATA
if (run$polls) {

  polls = read.csv("DataSets/polls.csv")
  dset.name = "polls"
  file.name = "bench-logitmm-polls.RData"

  idc  = !is.na(polls$bush)
  y    = polls$bush[idc]
  X.fe = model.matrix(bush ~ black + 0, data=polls[idc,])
  X.re = model.matrix(bush ~ state + 0, data=polls[idc,])
  n    = rep(1, length(y))
  
  P.a = ncol(X.re);
  P.b = ncol(X.fe);

  shape = 2
  rate  = 2
  
  prec.b = rep(0.01, P.b)
  
  m.0   = rep(0, P.b);
  P.0   = diag(prec.b, P.b);

  ## samp = 10000
  ## burn = 2000
  ## verbose = 1000
  ## df = Inf
  ## ntrials = 1

  bench.polls = list()
  if (load.old) load(file.name);
  
  for (i in run.methods) {
    nm = mm.method[i]
    bench.polls[[nm]] <- benchmark.blogit.mm(y, X.re, X.fe, n, shape, rate, m.0=m.0, C.0=solve(P.0),
                                             samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                             method=nm,
                                             dset.name=dset.name, df=df, var.names=var.names[[nm]], kappa=1)
  }

  polls.tbl = setup.table(bench.polls, "ab")

  if (write.it) save(y, X.re, X.fe, n, shape, rate, prec.b, m.0, P.0, df, bench.polls, polls.tbl, file=file.name);

  ## .........................................................
  ## Check sensitivity to hyperparameters and model structure.
  
  ## Without mean of random effects.
  check <- blogit.mm.map(y, X.re, X.fe, n, shape=shape, rate=rate, 
                         m.0=m.0, P.0=P.0, abphi.0=NULL, calc.V=TRUE, trace=FALSE)
  c(check$m[1:49], check$m[50])

  ## With mean of random effects.
  ## check.alt <- blogit.mm.map.alt(y, X.re, X.fe, n, shape=shape, rate=rate, kappa=0.1,
  ##                            m.0=m.0, P.0=P.0, abphim.0=NULL, calc.V=TRUE, trace=FALSE)
  ## c(check.alt$m[1:49], check.alt$m[50], check.alt$m[52])

  kappa = 100
  check.4 <- blogit.mm.map.4(y, X.re, X.fe, n, shape=shape, rate=rate, kappa=kappa,
                             m.0=m.0, P.0=P.0, abphim.0=NULL, calc.V=TRUE, trace=FALSE, maxit=10000)
  
  check.5 <- blogit.mm.map.5(y, X.re, X.fe, n, shape=shape, rate=rate, kappa=kappa,
                             m.0=m.0, P.0=P.0, abtm.0=NULL, calc.V=TRUE, trace=FALSE, maxit=10000)

  ## Check results using lme4 package.  I'm still unclear on their priors.
  hglm1 = glmer(bush ~ black + (1 | state), family=binomial(link="logit"), data=polls)

  coef(hglm1)
  r = ranef(hglm1, postVar = TRUE)
  summary(hglm1)
  ## Grab the posterior variances associated with each random effect
  attr(r$state,"postVar")
  
}

##------------------------------------------------------------------------------
## XEROP DATA
if (run$xerop) {
  
  require("epicalc")
  dset.name = "Xerop"
  file.name = "bench-logitmm-xerop.RData"
  
  data(Xerop)
  Xerop$id = factor(Xerop$id)

  y    = Xerop$respinfect
  X.fe = model.matrix(respinfect ~ xerop + age.month + sex + ht.for.age + stunted + factor(season) + 0, data=Xerop)
  X.re = model.matrix(respinfect ~ id + 0, data=Xerop)
  n    = rep(1, length(y))
  
  P.a = ncol(X.re);
  P.b = ncol(X.fe);

  shape = 2
  rate  = 2
  
  prec.b = rep(0.01, P.b)
  
  m.0   = rep(0, P.b);
  P.0   = diag(prec.b, P.b);
  n = rep(1, length(y))
  
  ## samp = 10000
  ## burn = 2000
  ## verbose = 100
  ## df = Inf
  ## ntrials = 1

  bench.xerop = list()
  if (load.old) load(file.name);
  
  for (i in run.methods) {
    ## Calculating Hessian numerically is a problem for a problem of this size.
    nm = mm.method[i]
    bench.xerop[[nm]] <- benchmark.blogit.mm(y, X.re, X.fe, n, shape, rate, m.0=m.0, C.0=solve(P.0),
                                             samp=samp, burn=burn, ntrials=ntrials, verbose=verbose,
                                             method=nm,
                                             dset.name=dset.name, df=df, var.names=var.names[[nm]], kappa=0)
  }

  xerop.tbl = setup.table(bench.xerop, "ab")

  if (write.it) save(y, X.re, X.fe, n, shape, rate, prec.b, m.0, P.0, df, bench.xerop, xerop.tbl, file=file.name);

  ## .........................................................
  hglm1 <- glmer(respinfect ~ xerop + age.month + sex + ht.for.age + stunted + factor(season) + (1 | id),
                 family=binomial(link="logit"), data=Xerop)
  ## This doesn't seem to work for glms.
  ## lme4.samp = mcmcsamp(hglm1, n = 12000)
  
}

if (FALSE) {

  ## FOR NEG BINOM

  springbok = read.csv("DataSets/springbok.csv", header=TRUE)

  ## Reset the baseline year
  springbok$timesince1990 = springbok$year - 1990
  ## Recast site as a factor
  springbok$site = factor(springbok$site)
  xyplot(counts ~ timesince1990 | site, data=springbok)
  hglm1 = glmer(counts ~ timesince1990 + (1+timesince1990 | site), data=springbok, family=poisson)

  y = springbok$counts
  
}

################################################################################
                                 ## APPENDIX ##
################################################################################

## I have encountered multiple difficulties that illuminate why inference for
## mixed models may be difficult.  First, I have found that there is sensitivity
## to the choice of prior and to subtle differences in the structure of the
## model.  In the mixed model synthetic data I created, when I use the second
## blogit routine, I find that the choice of shape, rate, and kappa plays an
## important role in the posterior mode.  Second, in the Xerop data set, when I
## try to calculate the hessian at the mode numerically I run into problems.
## There are like 283 dimensions, so maybe that is part of the problem.  This
## shows why it is really useful to use a sparse solver.  Most of the 283
## dimensional matrix will be zeros in the models we use above.  We could take
## advantage of this structure in the PG routine as well.  Third, it seems that
## the solution may be sensitive to the stopping conditions one uses.  Using the
## synthetic dataset I get noticably different answers for blogit.mm.llh.3 using
## factr=1e7 and factr=1e12 using the L-BFGS-B algorithm when kappa = 0.001.
## Another thing to remember is that transforming the prior changes the mode.
## Thus if one phrases things in terms of the log-precisions (log phi) instead
## of the precision (phi), then one arrives at different posterior modes.  (This
## is a reminder to me of why loss functions are so important.  Ultimately the
## loss function tells you about the coordinate system you are interested in.)

## It appears that I haven't been monitoring the convergence and that sometimes
## I have reached the max iter and stopped before convergence.  However, now
## that I have fixed that, even when I do get convergence I sometimes do not
## have identical answers.  One thing about transformation is that it will
## change how fast a value changes.  If I am working with phi, then phi changes
## more quickly when phi > 1 than if I was working with the same values on the
## log(phi) scale.  I suppose locally that is a linear relationship.
