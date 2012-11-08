
if (!is.loaded("BayesLogit.so")) dyn.load("BayesLogit.so")
else { dyn.unload("BayesLogit.so"); dyn.load("BayesLogit.so"); }
source("../R/LogitWrapper.R", local="TRUE");

out.D = rpg.devroye(100, 1.0, 0.0);

source("../R/PG.R");
source("../R/logit-MCMC.R");
source("../R/logit-EM.R");
source("../R/Multinomial.R");

# If you want to generate some data.
N = 1000;
x = cbind(1, rnorm(N));
## x = matrix(rnorm(2*N), N, 2)
beta = c(2.0, 1.0);
psi = x %*% as.matrix(beta);
p = 1 / (1 + exp(-1.0 * psi))
y = rbinom(N, 1, p);
n = rep(1, N);
y.prior = 0.5;
x.prior = colMeans(x);
n.prior = 1.0;

glm.1 = glm(y~x+0, family=binomial(link="logit"))
summary(glm.1)

## logit.EM.R(y,x,n)

out = logit(y, x, samp=1000, burn=100);
out.C = mlogit(y, x, samp=1000, burn=100);
out.R = mlogit.R(y, x, n = rep(1, length(y)), samp=100, burn=100);

## Multinomial data.
# If you want to generate some data.
N = 1000;
P = 4;
J = 4;
x = cbind(1, matrix(rnorm(N*(P-1)), N, P-1));
## x = matrix(rnorm(2*N), N, 2)
beta = matrix(rnorm(P*(J-1)), P, J-1);
psi = x %*% as.matrix(beta);
free.p = exp(psi) / (1 + rowSums(exp(psi)))
last.p = 1 - rowSums(free.p);
p = cbind(free.p, last.p);
y = array(0, dim=c(N, J));
for (i in 1:N)
  y[i,] = rmultinom(1, 1, p[i,]);
n = rep(1, N);
y = y[,1:(J-1)];

start = proc.time();
out.R = mlogit.R(y, x, samp=1000, burn=100, verbose=1000);
end = proc.time(); print(end-start);
start = proc.time();
out.C = mlogit(y, x, samp=1000, burn=100);
end = proc.time(); print(end-start);

## # What the original R example uses.
## x = c(53,56,57,63,66,67,68,69, 70,72,73, 75,76,78,79,80,81)
## y = c( 1, 1, 1, 0, 0, 0, 0, 0,3/4, 0, 0,1/2, 0, 0, 0, 0, 0)
## n = c( 1, 1, 1, 1, 1, 3, 1, 1,  4, 1, 1,  2, 2, 1, 1, 1, 1)

## y.prior = 0;
## x.prior = 81;
## n.prior = 1;

data("spambase")
y = spambase$is.spam;
x = spambase[c("word.freq.free", "char.freq.!", "capital.run.length.average")];

start = proc.time();
# out.R = logit.R(y, x, samp=1000, burn=100);
out.R = logit.R(c(y.prior,y), rbind(x.prior, x), c(n.prior, n), 1000, 500);
R.time = proc.time() - start;
print(R.time);

#out.C = logit(y, x, n, y.prior, x.prior, n.prior, 10000, 500);
out.C = logit(y, x, samp=1000, burn=500);

out.glm = glm(y ~ as.matrix(x)-1, family=binomial(link="logit"));

print(colMeans(out.R$beta));
print(colMeans(out.C$beta));

t1 = proc.time();
out.R = rpg.R(10000, 1.0, 1.0);
t2 = proc.time();
print(t2-t1);
# out.P = c(); for(i in 1:10000) { out.P[i] = rpg.P (0.0); }
t1 = proc.time();
out.C = rpg  (100000, 1.0, 1.0);
t2 = proc.time();
print(t2-t1)
t1 = proc.time()
out.D = rpg.devroye(100000, 1, 1.0);
t2 = proc.time()
print(t2-t1)

print(c(mean(out.R), sd(out.R), mean(out.R^2), mean(out.R^3), mean(out.R^4)));
print(c(mean(out.C), sd(out.C), mean(out.C^2), mean(out.C^3), mean(out.C^4)));
print(c(mean(out.D), sd(out.D), mean(out.D^2), mean(out.D^3), mean(out.D^4)));
# print(c(mean(out.P), sd(out.P), mean(out.P^2), mean(out.P^3), mean(out.P^4)));

################################################################################
                            ## DATA CREATION HELP ##
################################################################################

# This is the Wisconsin Breast Cancer data.
wdbc = read.csv("wdbc.data.txt", header=FALSE);
# wpbc = read.csv("wpbc.data.txt", header=FALSE);

features=c("radius",
           "texture",
           "perimeter",
           "area",
           "smoothness",
           "compactness",
           "concavity",
           "concave points",
           "symmetry",
           "fractal.dimension");

all.measures = c();
for(i in 1:10){
    all.measures = c(all.measures,
                     paste(features[i], ".mean", sep=""),
                     paste(features[i], ".se", sep=""),
                     paste(features[i], ".worst", sep=""));
}

the.names = c("ID", "Diagnosis", all.measures);
names(wdbc) = the.names;
# names(wpbc) = the.names;

# Change "M" and "B" to 1 and 0 respectively.
wdbc$Diagnosis = (wdbc$Diagnosis=="M")*1

## # We need to standardize x, otherwise x'x is almost singular.
## for(i in 3:32){
##     wdbc[,i] = (wdbc[,i] - mean(wdbc[,i]));
##     wdbc[,i] = wdbc[,i] / sqrt(sum(wdbc[,i]^2));
## }

save(wdbc, file="wdbc.RData");

## only.the.means = c(1,2,seq(3, 32, 3));
## wdbc = wdbc[,only.the.means];
## wpbc = wpbc[,only.the.means];

# The spambase database.
spam = read.csv("spambase.data.txt", header=FALSE);
the.names = scan("spambase.names", what="character");
names(spam) = the.names;

spambase = cbind("is.spam"=spam$spam, spam[,1:57]);

save(spambase, file="spambase.RData");

################################################################################

## SCRATCH

