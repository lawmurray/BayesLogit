#!/bin/bash

## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

FILENAME=$0
RELDIR=${FILENAME%MakeDynSupplement.bash}
BASE=${RELDIR} # Code/R directory.
DEST="$BASE/../DynSupplement"
DATASOURCE="Benchmark-DataSets"
DATADEST="$DEST/$DATASOURCE"

#echo "File name is $FILENAME"
#echo "Relative directory is $RELDIR"
#echo "base directory is $BASE"

rsyncit="rsync -Crzut $1 --exclude-from=$HOME/.rsync-exclude"

# R Files
RFILES=( 
	"Benchmark-Utilities.R"
	"Benchmark-DynNB.R"
	"DynNBPG.R"
	"DynNBFS-2009.R"
	"DynNBCUBS.R"
	"DynNBOmegaBlock.R"
	"Benchmark-DynLogit.R"
	"DynLogitPG.R"
	"DynLogitdRUM.R"       
	"DynLogitCUBS.R"
	"DynLogitOmegaBlock.R"
	"Stationary.R"
	"compmix.R"
	"DynLogitBetaBlockMH.R"
	"NB-Shape.R"
)

# Data Files
DFILES=(
	"DynLogit-synth-high-2-n-1.RData"
	"DynLogit-synth-high-2-n-20.RData"
	"DynLogit-synth-low-2-n-1.RData"
	"DynLogit-synth-low-2-n-20.RData"
	"DynNB-synth-high-2-mu-100.RData"
	"DynNB-synth-high-2-mu-10.RData"
	"DynNB-synth-low-2-mu-100.RData"
	"DynNB-synth-low-2-mu-10.RData"
)
 
# Print out array.
# echo ${RFILES[*]}

for file in ${RFILES[*]}
do
	# echo ${BASE}$file
	$rsyncit ${BASE}$file $DEST/
done

for file in ${DFILES[*]}
do
	# echo ${BASE}$DATASOURCE/$file
	$rsyncit ${BASE}$DATASOURCE/$file $DATADEST/
done