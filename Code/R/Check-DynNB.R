## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


data("rain", package="BayesLogit")
load("bench-dynlogit-tokyo.RData")

################################################################################
                             ## CHECK SYNTHETIC ##
################################################################################

P = 2
nb.mean = 10
corr.type = "high"
## est.ar = "with.ar"
est.ar = "wout.ar"

## for (est.ar in c("wout.ar", "with.ar")) {
## for (P in c(2,4)) {
## for (nb.mean in c(10,100)) {
## for (corr.type in c("low", "high")) {

dset.name = paste(corr.type, "-", P, "-mu-", nb.mean, sep="");
source.file = paste("DynNB-synth-", dset.name, ".RData", sep="")
load(file.path("Benchmark-DataSets", source.file))
filename = paste("bench-dynnb-", dset.name, "-", est.ar, ".RData", sep="")
load(filename)

psi.pg = mean(bench.synth$PG$gb$alpha) + colSums(synth.table$ave.sstat[,-1,1,1] * t(X))
psi.fs = mean(bench.synth$FS$gb$alpha) + colSums(synth.table$ave.sstat[,-1,1,2] * t(X))
psi.cubs = mean(bench.synth$CUBS$gb$alpha) + colSums(synth.table$ave.sstat[,-1,1,3] * t(X))

## psi.pg = colSums(synth.table$ave.sstat[,-1,1,1] * t(X))
## psi.fs = colSums(synth.table$ave.sstat[,-1,1,2] * t(X))

par(mfrow=c(3,1))

plot(beta[1,], type="l");
lines(synth.table$ave.sstat[1,,1,1], col=2)
lines(synth.table$ave.sstat[1,,1,2], col=3)
lines(synth.table$ave.sstat[1,,1,3], col=4)
lines(bench.synth$PG$gb$beta[1000,1,], col=5)

plot(beta[2,], type="l");
lines(synth.table$ave.sstat[2,,1,1], col=2)
lines(synth.table$ave.sstat[2,,1,2], col=3)
lines(synth.table$ave.sstat[2,,1,3], col=4) 

plot(log.mean, type="l")
lines(psi.pg, col=2)
lines(psi.fs, col=3)
lines(psi.cubs, col=4)

par(mfrow=c(1,1))
for (i in seq(1000, 5000, 100)) {
  plot(beta[1,], type="l");
  lines(synth.table$ave.sstat[1,,1,1], col=rgb(1, 0, 0, 0.2))
  lines(synth.table$ave.sstat[1,,1,2], col=rgb(0, 1, 0, 0.2))
  lines(synth.table$ave.sstat[1,,1,3], col=rgb(0, 0, 1, 0.2))
  lines(bench.synth$PG$gb$beta[i,1,], col=rgb(1, 0, 0, 1))
  lines(bench.synth$FS$gb$beta[i,1,], col=rgb(0, 1, 0, 1))
  lines(bench.synth$CUBS$gb$beta[i,1,], col=rgb(0, 0, 1, 1))
  readline("<ENTER>")
}
