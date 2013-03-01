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

set +x

FILES="CUBS.Rd AR1.Rd AR1Ind.Rd"

## Add files for dynamic.
if [[ $1 == Dynamic ]]; then
  for FILE in $FILES; do
	set -x
    cp "${RELDIR}${1}/$FILE" "${RELDIR}../man/$FILE"
	set +x
  done
fi

## Remove files for static.
if [[ $1 == Static ]]; then
  for FILE in $FILES; do
	PTF="${RELDIR}../man/$FILE"
	if [[ -e "$PTF" ]]; then
	  set -x
	  rm "$PTF"
	  set +x
    fi
  done	
fi
