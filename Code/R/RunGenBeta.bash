## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

## Number of files
rd=$(ls *.RData)
rd=( $rd )
N=${#rd[@]}

echo "Number of files: $N"
let "i = N+1"

while ((N < 10)); do
    let "N=$N+1"
    echo "Working on ${N}."
    R CMD BATCH --no-save GenerateBeta.R
    mv GenerateBeta.Rout "GenerateBeta.Rout.$i"
    rd=$(ls *.RData)
    rd=( $rd )
    N=${#rd[@]}
done