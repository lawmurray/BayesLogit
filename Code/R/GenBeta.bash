## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

filename="GenerateBeta${1}.R"
clerror="0"

## Value and number of arguments.  Starts at 0, as opposed to $1, $2, etc.
ARGV=( $@ )
ARGC=${#ARGV[@]}

## Check for command line errors.

## Exit if file does not exist.
if ! [[ -e "$filename" ]] ; then
    echo "Could not find file <$filename>."
    clerror=1
fi

for((i=1;i<$ARGC;i++))
do
    if ! [[ "${ARGV[$i]}" =~ ^[0-9]*$ ]] ; then
	let "j=i+1"
	echo "Argument ${ARGV[$i]}: argument $j must be a number or nothing."
	clerror=1
    fi
done

## Exit if error.
if ! [[ $clerror -eq 0 ]]; then
    exit 15
fi

## Now run script.

echo "Will run <$filename>."

# Set available methods.
meth=(FS PG)
if [[ "$1" == "Logit" ]]; then
    meth=(FS PG HH DFAE) # For Logit
fi

# Number of methods.
let "N=${#meth[@]}-1"

tempfile="GenerateBetaTemp.R"
cp $filename $tempfile

## Change the samples, burn, or trials.
if ! [[ -z "$2" ]] ; then
    echo "Changing number of samples to $2."
    sed -i~ -e s/^samp\ \=\ [0-9]*/samp\ \=\ $2/ $tempfile
fi

if ! [[ -z "$3" ]] ; then
    echo "Changing number of burn-in to $3."
    sed -i~ -e s/^burn\ \=\ [0-9]*/burn\ \=\ $3/ $tempfile
fi

if ! [[ -z "$4" ]] ; then
    echo "Changing number of trials to $4."
    sed -i~ -e s/^trials\ =\ [0-9]*/trials\ =\ $4/ $tempfile
fi

# Go.
for i in $(seq 0 $N)
do
    M=${meth[i]}
    echo "Working on $M."
    ## sed -i~ -e s/FILL_IN/$M/ $tempfile
    sed -i~ -e s/^method\ =\ \".*\"/method\ =\ \"$M\"/ $tempfile

    R CMD BATCH --no-save $tempfile
    mv ${tempfile}out "${filename}out.$M"

done

##########

## Old Stuff ##

    # # One at a time.
    # ## Number of files
    # rd=$(ls *$M*.RData)
    # rd=( $rd )
    # N=${#rd[@]}

    # echo "Number of files: $N"
    # let "i = N+1"

    # while ((N < 10)); do
    # 	let "N=$N+1"
    # 	echo "Working on ${N}."
    # 	R CMD BATCH --no-save GenerateBetaTemp.R
    # 	mv GenerateBetaTemp.Rout "GenerateBeta.Rout.$M.$i"
    # 	rd=$(ls *$M*.RData)
    # 	rd=( $rd )
    # 	N=${#rd[@]}
    # done