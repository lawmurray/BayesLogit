FILENAME=$0

RELDIR=${FILENAME%set.bash}

## rsyncit="rsync -Crzut"

if [[ -z $1 ]]; then
    echo "You must specify \"Dynamic\" or \"Static\"."
    exit
fi

set -x

cp ${RELDIR}${1}/DESCRIPTION ${RELDIR}../DESCRIPTION
cp ${RELDIR}${1}/NAMESPACE   ${RELDIR}../NAMESPACE
cp ${RELDIR}${1}/Makevars    ${RELDIR}../src/Makevars
