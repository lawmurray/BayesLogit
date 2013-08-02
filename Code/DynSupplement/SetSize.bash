#!/bin/bash

# Default
samp=10000
burn=2000
verb=1000
reps=10
idcs="1:3"

# small
if [[ $1 == "small" ]]; then
	samp=100
	burn=20
	verb=10
	reps=2
fi

# large
if [[ $1 == "large" ]]; then
	samp=10000
	burn=2000
	verb=1000
	reps=10
fi

echo "Setting..."
echo "Number of samples: $samp"
echo "Burn in: $burn"
echo "Report mod: $verb"
echo "Batches: $reps"
echo "Will run all data sets."

sedexp='s/^samp = [0-9]*$/samp = '${samp}'/ 
        s/^burn = [0-9]*$/burn = '${burn}'/
        s/^verbose = [0-9]*$/verbose = '${verb}'/
        s/^ntrials = [0-9]* #/ntrials = '${reps}' #/
        s/"allsynth"=FALSE/"allsynth"=TRUE/
        s/run\.idc = .*$/run\.idc = '${idcs}'/'

# echo "$sedexp"
sed -i~ -e "$sedexp" Benchmark-DynLogit.R
sed -i~ -e "$sedexp" Benchmark-DynNB.R