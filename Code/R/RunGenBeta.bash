
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