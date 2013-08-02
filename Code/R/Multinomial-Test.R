## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


library("BayesLogit", lib.loc="~/RPackage/BayesLogit/Code/BLPackage/Test");
source("Multinomial.R")

## IRIS ##
##------------------------------------------------------------------------------
if (FALSE) {

  ## iris = read.csv("Datasets/iris.csv", header=TRUE)

  N = nrow(iris)
  P = ncol(iris)
  J = nlevels(iris$Species)

  X     = model.matrix(Species ~ ., data=iris);
  y.all = model.matrix(~ Species - 1, data=iris);
  y     = y.all[,-J];

  out = mlogit(y, X, rep(1, N), samp=1000);

  beta.ave = cbind(apply(out$beta, c(2,3), mean), 0);
  psi.ave  = X %*% beta.ave
  fit = apply(psi.ave, 1, which.max)
  ## View(fit)
  fit
  num.wrong = sum( abs(as.numeric(iris$Species) - fit) )

}

## WINE ##
##------------------------------------------------------------------------------
if (FALSE) {

  white = read.csv("Datasets/winequality-white.csv", sep=";");
  white$quality = factor(white$quality);
  
  N = nrow(white)
  P = ncol(white)
  J = nlevels(white$quality)

  X = as.matrix(white[,1:(P-1)]);
  y = matrix(0, N, J-1);

  for (i in 1:(J-1)) {
    y[,i] = white$quality==levels(white$quality)[i];
  }

  out = mult.logit.gibbs(y, X, rep(1, N), samp=1000);

  beta.ave = cbind(apply(out$beta, c(2,3), mean), 0);
  psi.ave  = X %*% beta.ave
  fit = apply(psi.ave, 1, which.max)
  ## View(fit)
  fit
  num.wrong = sum( abs(as.numeric(white$quality) - fit) )
  
}

## Synthetic
##------------------------------------------------------------------------------

if (FALSE) {

  ## This type of synthetic data may not be the best.  You may not predictors that are normal.
  
  N = 1000;
  P = 3;
  J = 2;
  num = 1;

  samp = 1000;

  beta = matrix(rnorm(P*(J-1)), P, J-1);
  X    = cbind(1, matrix(rnorm(N*(P-1), c(-2,-1,0,1,2), 1), N, P-1));
  psi  = X %*% beta;

  free.p = exp(psi) / (1 + rowSums(exp(psi)))
  last.p = 1 - rowSums(free.p);
  p      = cbind(free.p, last.p);

  ## True Multinomial Data.
  y.all = array(0, dim=c(N,J));
  for (i in 1:N)
    y.all[i,] = rmultinom(1, num, p[i,]);
  y = y.all[,-J] / N;
  n = rep(num, N);

  ## GLM
  ## glm.1 = glm(response ~ predictors, family=binomial, data=a.df);
  ## summary(glm.1)
  
  ## BayesLogit
  start = proc.time()
  out.B = mlogit(y, X, samp=samp)
  end = proc.time()

  print(end-start);

}

## Synthetic
##------------------------------------------------------------------------------
if (FALSE) {

  N = 1000;
  P = 4;
  J = 3;
  num = 1;

  samp = 1000;

  beta = matrix(rnorm(P*(J-1)), P, J-1);
  X    = cbind(1, matrix(
                  rnorm(N*(P-1), c(-2,-1,0,1,2), 1),
                  N, P-1));
  psi  = X %*% beta;

  free.p = exp(psi) / (1 + rowSums(exp(psi)))
  last.p = 1 - rowSums(free.p);
  p      = cbind(free.p, last.p);

  ## Only one roll per ith sample.
  y.cat = rep(0, N);
  for (i in 1:N)
    y.cat[i] = sample.int(J, 1, prob=p[i,]);
  y.all = model.matrix( ~ factor(y.cat) - 1 );
  y = y.all[,-J];
  n = rep(1, N);

  a.df = data.frame("response"=y.cat, "predictors"=X[,-1]);

  ## ## True Multinomial Data.
  ## y.all = array(0, dim=c(N,J));
  ## for (i in 1:N)
  ##   y.all[i,] = rmultinom(num, 1, p[i,]);
  ## y = y.all[,-J]
  ## n = rep(num, N);

  ## GLM
  ## glm.1 = glm(response ~ predictors, family=binomial, data=a.df);
  ## summary(glm.1)
  
  ## BayesLogit
  start = proc.time()
  out.B = mlogit(y, X, samp=samp)
  end = proc.time()

  print(end-start);

  ## MCMCpack
  start = proc.time()
  ## mcmc.method=c("IndMH", "RWM", "slice")
  out.M = MCMCmnl(response ~ ., mcmc=samp, data=a.df, baseline=J, mcmc.method="slice")
  ## out.M = MCMCmnl(factor(y.cat) ~ X, mcmc=samp) 
  end = proc.time()

  print(end-start);

  ## Summary Statistics
  beta.mat = array(out.B$beta, dim=c(samp, P*(J-1)));
  mean(apply(beta.mat, 2, effectiveSize));
  mean(apply(beta.mat, 2, ESS));

  out.M.cube = array(0, dim=c(samp, P, J-1));
  for (i in 1:samp) out.M.cube[i,,] = t( array(out.M[i,], dim=c(J-1, P)) )
  
  mean(effectiveSize(out.M))
  mean(apply(out.M, 2, ESS));

  beta.B = apply(out.B$beta, c(2,3), mean)
  beta.M = apply(out.M.cube, c(2,3), mean)

  psi.B = X %*% beta.B
  psi.M = X %*% beta.M
  p.B = exp(psi.B) / (1 + rowSums(exp(psi.B)))
  p.M = exp(psi.M) / (1 + rowSums(exp(psi.M)))

  ## Check by visual
  for (j in 1:(J-1)) {
    for (i in 1:P) {
      par(mfrow=c(1,2));
      hist(out.B$beta[,i,j] - beta[i,j], breaks=20)
      hist(out.M.cube[,i,j] - beta[i,j], breaks=20)
      cat(sd(out.B$beta[,i,j]), sd(out.M.cube[,i,j]), " ");
      readline("PRESS ENTER");
    }
  }
  
}

## NETHVOTE ##
##------------------------------------------------------------------------------
if (FALSE) {

  source("~/RPackage/BayesLogit/Code/R/Efficiency.R")
  
  library("BayesLogit")
  library("MCMCpack");

  data("Nethvote");
  samp = 1000;
  
  form = vote ~ . -distPvdA

  the.levels = levels(factor(Nethvote$vote));
  J = length(the.levels);
  
  ## samp = 1000;
  ## K = 100;
  ## ss = sample.int(N, K);

  nv = Nethvote
  ## nv = Nethvote[ss,]
  nv$vote = factor(nv$vote);
  
  y.all = model.matrix(~vote-1, data=nv)
  y = y.all[,1:(J-1)]
  X = model.matrix(form, data=nv)
  N = nrow(X);
  P = ncol(X);
  n = rep(1, N);
 
  start = proc.time()
  out.B = mlogit(y, X, samp=samp)
  end = proc.time()

  print(end-start);

  start = proc.time()
  out.M = MCMCmnl(form, data=nv, mcmc.method="IndMH", baseline=the.levels[J], mcmc=samp) 
  end = proc.time()

  print(end-start);

  beta.mat = array(out.B$beta, dim=c(samp, 30));
  mean(apply(beta.mat, 2, effectiveSize));
  mean(apply(beta.mat, 2, ESS));
  mean(effectiveSize(out.M))
  mean(apply(out.M, 2, ESS));

  out.M.cube = array(0, dim=c(samp, P, J-1));
  for (i in 1:samp) out.M.cube[i,,] = t( array(out.M[i,], dim=c(J-1, P)) )

  beta.B = apply(out.B$beta, c(2,3), mean) 
  beta.M = colMeans(out.M.cube);
  
}

## Glass ##
##------------------------------------------------------------------------------
if (FALSE) {

  source("~/RPackage/BayesLogit/Code/R/Efficiency.R")
  
  library("BayesLogit")
  library("MCMCpack");

  glass = read.csv("~/RPackage/BayesLogit/Code/R/DataSets/glass.dat", header=FALSE)
  colnames(glass) = c("Id", "RI", "Na", "Mg", "Al", "Si", "K", "Ca", "Ba", "Fe", "class");
  glass = glass[,2:11];

  form = class ~ RI + Na + Mg + Al + Si + K + Ca + Ba + Fe
  ## form = class ~ RI + Al + Si + K + Ca + Fe
  ## form = class ~ RI + Al + Si + K + Ca
  ## form = class ~ RI

  the.levels = levels(factor(glass$class));
  J = length(the.levels);

  samp = 10000;
  burn = 1000;
  
  glass$class = factor(glass$class);
  
  y.all = model.matrix(~class-1, data=glass)
  y = y.all[,1:(J-1)]
  X = model.matrix(form, data=glass)
  N = nrow(X);
  P = ncol(X);
  n = rep(1, N)

  p.0 = 0;
  
  ## BayesLogit
  P.0 = array(diag(p.0, P), dim=c(P, P, J-1));
  start = proc.time()
  out.B = mlogit(y, X, samp=samp, burn=burn, P.0=P.0)
  end = proc.time()

  print(end-start);

  ## MCMCpack
  start = proc.time()
  out.M = MCMCmnl(form,
    data=glass, mcmc=samp, mcmc.method="IndMH", baseline=the.levels[J], B0=p.0) 
  end = proc.time()

  print(end-start);

  ## Organize Output
  beta.mat = array(out.B$beta, dim=c(samp, P*(J-1)));
  mean(apply(beta.mat, 2, effectiveSize));
  mean(apply(beta.mat, 2, ESS));
  mean(effectiveSize(out.M))
  mean(apply(out.M, 2, ESS));

  out.M.cube = array(0, dim=c(samp, P, J-1));
  for (i in 1:samp) out.M.cube[i,,] = t( array(out.M[i,], dim=c(J-1, P)) )

  beta.B = apply(out.B$beta, c(2,3), mean)
  beta.M = apply(out.M.cube, c(2,3), mean)

  psi.B = X %*% beta.B
  psi.M = X %*% beta.M
  p.B = exp(psi.B) / (1 + rowSums(exp(psi.B)))
  p.B = cbind(p.B, 1 - rowSums(p.B))
  p.M = exp(psi.M) / (1 + rowSums(exp(psi.M)))

  ## Let's see how many we got right.
  guess = apply(p.B, 1, which.max);
  correct = as.numeric(glass$class==the.levels[guess]);
  compare = cbind(glass$class, guess, correct);

  rate = array(0, dim=c(J,3));

  for (i in 1:J) {
    cl  = glass$class==the.levels[i]
    cor = correct[cl]
    rate[i,] = c(sum(cor), sum(cl), sum(cor) / sum(cl))
  }

  ## Log-odds
  psi.samp = array(0, dim=c(samp, N, J-1));
  for (i in 1:samp) {
    psi.samp[i,,] = X %*% out.B$beta[i,,];
  }

  Q = quantile(psi.samp[,1,1], c(0.05, 0.95));
  plot( rep(1+.5/J, 2), Q, type="l", xlim=c(1,N), ylim=c(-200, 100));
  for (i in 1:N) {
    for (j in 1:(J-1)) {
      Q = quantile(psi.samp[,i,j], c(0.05, 0.95));
      lines( rep(i+.5*j/J, 2), Q, col=j);
    }
  }

  
  
}
