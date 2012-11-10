#!/bin/bash
# Copy the necessary files to the BayesLogit Package directory.

## BASE=THE DIRECTORY WHERE YOU GIT CLONED BAYESLOGIT.
BASE=/Users/jwindle/RPackage/BayesLogit

rsyncit="rsync -Crvzut --exclude-from=$BASE/.rsync-exclude"

CODE=$BASE/Code
INCL=$CODE/C/include

BLDIR=$CODE/BLPackage/BayesLogit

# CPP files.
$rsyncit $CODE/C/Logit.hpp        $BLDIR/src/Logit.h
$rsyncit $CODE/C/MultLogit.hpp    $BLDIR/src/MultLogit.h
$rsyncit $CODE/C/LogitWrapper.hpp $BLDIR/src/LogitWrapper.h
$rsyncit $CODE/C/LogitWrapper.cpp $BLDIR/src/
$rsyncit $CODE/C/PolyaGamma.hpp   $BLDIR/src/PolyaGamma.h
$rsyncit $CODE/C/FFBS.h           $BLDIR/src/
$rsyncit $CODE/C/FFBS.cpp         $BLDIR/src/
$rsyncit $CODE/C/FSF_nmix.hpp     $BLDIR/src/FSF_nmix.h
$rsyncit $CODE/C/FSF_nmix.cpp     $BLDIR/src/

$rsyncit $INCL/Matrix/Matrix.h             $BLDIR/src/
$rsyncit $INCL/Matrix/Matrix.cpp           $BLDIR/src/
$rsyncit $INCL/Matrix/MatrixFrame.h        $BLDIR/src/
$rsyncit $INCL/Matrix/MatrixFrame.cpp      $BLDIR/src/
$rsyncit $INCL/RNG/RNG.hpp                 $BLDIR/src/RNG.h
$rsyncit $INCL/RNG/RNG.cpp                 $BLDIR/src/
$rsyncit $INCL/RNG/RRNG.hpp                $BLDIR/src/RRNG.h
$rsyncit $INCL/RNG/RRNG.cpp                $BLDIR/src/
$rsyncit $INCL/Normal.hpp                  $BLDIR/src/Normal.h

# R files.
$rsyncit $CODE/R/LogitWrapper.R     $BLDIR/R/
$rsyncit $CODE/R/NB-Indicators.R    $BLDIR/R/
$rsyncit $CODE/R/Logit-Indicators.R $BLDIR/R/
$rsyncit $CODE/R/FFBS.R             $BLDIR/R/
$rsyncit $CODE/R/KS.R               $BLDIR/R/
## The function calls in these files are not in the NAMESPACE.
$rsyncit $CODE/R/LogitPG.R       $BLDIR/R/
$rsyncit $CODE/R/logit-EM.R      $BLDIR/R/
$rsyncit $CODE/R/logit-combine.R $BLDIR/R/
$rsyncit $CODE/R/PG.R            $BLDIR/R/
$rsyncit $CODE/R/MultLogitPG.R   $BLDIR/R/

# Data files.
$rsyncit $CODE/R/DataSets/spambase.RData      $BLDIR/data/

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