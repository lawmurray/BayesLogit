FILENAME=$0

RELDIR=${FILENAME%set.bash}

## rsyncit="rsync -Crzut"

if [[ -z $1 ]]; then
    echo "You must specify \"Dynamic\" or \"Static\"."
    exit
fi

name=`echo $1 | tr '[:upper:]' '[:lower:]'` 
name=`echo $name | sed -e 's/^./\U&/'`
#echo $name

set -x

cp ${RELDIR}${name}/DESCRIPTION ${RELDIR}../DESCRIPTION
cp ${RELDIR}${name}/NAMESPACE   ${RELDIR}../NAMESPACE
cp ${RELDIR}${name}/Makevars    ${RELDIR}../src/Makevars

set +x

FILES="man/CUBS.Rd man/AR1.Rd man/AR1Ind.Rd R/CUBS.R"

## Add files for dynamic.
if [[ $name == Dynamic ]]; then
  for FILE in $FILES; do
	set -x
    cp "${RELDIR}${name}/$FILE" "${RELDIR}../$FILE"
	set +x
  done
fi

## Remove files for static.
if [[ $name == Static ]]; then
  for FILE in $FILES; do
	PTF="${RELDIR}../$FILE"
	if [[ -e "$PTF" ]]; then
	  set -x
	  rm "$PTF"
	  set +x
    fi
  done	
fi
