#!/bin/bash
# Copy the necessary files to the BayesLogit Package directory.

rsyncit="rsync -Crvzut --exclude-from=$HOME/.rsync-exclude"

BASE=/Users/jwindle/RPackage/BayesLogit/Code
MYLIB=/Users/jwindle/RV-Project/Code/C_Examples/MyLib

BLDIR=/Users/jwindle/RPackage/BayesLogit/Code/BLPackage/BayesLogit

# CPP files.
$rsyncit $BASE/C/Logit.hpp        $BLDIR/src/Logit.h
$rsyncit $BASE/C/MultLogit.hpp    $BLDIR/src/MultLogit.h
$rsyncit $BASE/C/LogitWrapper.hpp $BLDIR/src/LogitWrapper.h
$rsyncit $BASE/C/LogitWrapper.cpp $BLDIR/src/
$rsyncit $BASE/C/PolyaGamma.hpp   $BLDIR/src/PolyaGamma.h

$rsyncit $MYLIB/Matrix/Matrix.h             $BLDIR/src/
$rsyncit $MYLIB/Matrix/MatrixFrame.h        $BLDIR/src/
$rsyncit $MYLIB/RNG/RNG.hpp                 $BLDIR/src/RNG.h
$rsyncit $MYLIB/RNG/RRNG.hpp                $BLDIR/src/RRNG.h
$rsyncit $MYLIB/RandomVariates/Normal.hpp   $BLDIR/src/Normal.h

# R files.
$rsyncit $BASE/C/LogitWrapper.R  $BLDIR/R/
$rsyncit $BASE/R/KS.R            $BLDIR/R/
## The function calls in these files are not in the NAMESPACE.
$rsyncit $BASE/R/LogitPG.R    $BLDIR/R/
$rsyncit $BASE/R/logit-EM.R      $BLDIR/R/
$rsyncit $BASE/R/logit-combine.R $BLDIR/R/
$rsyncit $BASE/R/PG.R            $BLDIR/R/
$rsyncit $BASE/R/MultLogitPG.R   $BLDIR/R/

# Data files.
$rsyncit $BASE/R/DataSets/spambase.RData      $BLDIR/data/

sed -i~ s/\"Logit\"/\"BayesLogit\"/ $BLDIR/R/LogitWrapper.R

# Change to .h
sed -i~ s/\.hpp/.h/g $BLDIR/src/*.[ch]
sed -i~ s/\.hpp/.h/g $BLDIR/src/*.cpp

# Change to Rprintf
sed -i~ s/fprintf\(stderr,/printf\(/g $BLDIR/src/*.[ch]
sed -i~ s/fprintf\(stderr,/printf\(/g $BLDIR/src/*.cpp

# There must be a better way.
sed -i~ -e 's/Rprintf/printf/g' $BLDIR/src/*.[ch]
sed -i~ -e 's/Rprintf/printf/g' $BLDIR/src/*.cpp

sed -i~ -e 's/printf/Rprintf/g' $BLDIR/src/*.[ch]
sed -i~ -e 's/printf/Rprintf/g' $BLDIR/src/*.cpp