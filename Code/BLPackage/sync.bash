#!/bin/bash
# Copy the necessary files to the BayesLogit Package directory.
# use -v flag, e.g. bash sync.bash -v for verbose.

# Find the directory sync.bash is in, then find base directory.
FILENAME=$0
RELDIR=${FILENAME%sync.bash}
BASE=${RELDIR}../..
# echo "File name is $FILENAME"
# echo "Relative directory is $RELDIR"
# echo "base directory is $BASE"

rsyncit="rsync -Crzut $1 --exclude-from=$BASE/.rsync-exclude"

CODE=$BASE/Code
INCL=$CODE/C/include

BLDIR=$CODE/BLPackage/BayesLogit

# CPP files.
$rsyncit $CODE/C/Logit.hpp         $BLDIR/src/Logit.h
$rsyncit $CODE/C/MultLogit.hpp     $BLDIR/src/MultLogit.h
$rsyncit $CODE/C/LogitWrapper.h    $BLDIR/src/
$rsyncit $CODE/C/LogitWrapper.cpp  $BLDIR/src/
$rsyncit $CODE/C/PolyaGamma.h      $BLDIR/src/
$rsyncit $CODE/C/PolyaGammaAlt.h   $BLDIR/src/
$rsyncit $CODE/C/PolyaGammaSP.h    $BLDIR/src/
$rsyncit $CODE/C/PolyaGamma.cpp    $BLDIR/src/
$rsyncit $CODE/C/PolyaGammaAlt.cpp $BLDIR/src/
$rsyncit $CODE/C/PolyaGammaSP.cpp  $BLDIR/src/
$rsyncit $CODE/C/InvertY.hpp       $BLDIR/src/InvertY.h
$rsyncit $CODE/C/InvertY2.hpp       $BLDIR/src/InvertY2.h
$rsyncit $CODE/C/InvertY.cpp       $BLDIR/src/
$rsyncit $CODE/C/InvertY2.cpp       $BLDIR/src/
## $rsyncit $CODE/C/FFBS.h           $BLDIR/src/
## $rsyncit $CODE/C/FFBS.cpp         $BLDIR/src/
$rsyncit $CODE/C/FSF_nmix.hpp     $BLDIR/src/FSF_nmix.h
$rsyncit $CODE/C/FSF_nmix.cpp     $BLDIR/src/
$rsyncit $CODE/C/CUBS.h           $BLDIR/src/
$rsyncit $CODE/C/CUBS.cpp         $BLDIR/src/
$rsyncit $CODE/C/CUBS_update.h    $BLDIR/src/
$rsyncit $CODE/C/CUBS_update.cpp  $BLDIR/src/
$rsyncit $CODE/C/AR1.h            $BLDIR/src/
$rsyncit $CODE/C/AR1.cpp          $BLDIR/src/
$rsyncit $CODE/C/DynExpFamMH.hpp  $BLDIR/src/DynExpFamMH.h
$rsyncit $CODE/C/DynExpFamMH.cpp  $BLDIR/src/
  

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
$rsyncit $CODE/R/Indicators.R       $BLDIR/R/
$rsyncit $CODE/R/Logit-Indicators.R $BLDIR/R/
$rsyncit $CODE/R/FFBS.R             $BLDIR/R/
$rsyncit $CODE/R/CUBS.R             $BLDIR/inst/Dynamic/R/
$rsyncit $CODE/R/DynExpFamMHWrapper.R $BLDIR/R/
$rsyncit $CODE/R/KS.R               $BLDIR/R/
$rsyncit $CODE/R/AR1.R              $BLDIR/R/
$rsyncit $CODE/R/ComputeMixture.R   $BLDIR/R/
$rsyncit $CODE/R/compmix.R          $BLDIR/R/
## The function calls in these files are not in the NAMESPACE.
$rsyncit $CODE/R/LogitPG.R       $BLDIR/R/
$rsyncit $CODE/R/logit-EM.R      $BLDIR/R/
$rsyncit $CODE/R/logit-combine.R $BLDIR/R/
$rsyncit $CODE/R/PG.R            $BLDIR/R/
$rsyncit $CODE/R/MultLogitPG.R   $BLDIR/R/


# Data files.
$rsyncit $CODE/R/DataSets/spambase.RData   $BLDIR/data/
$rsyncit $CODE/R/DataSets/rain.RData       $BLDIR/data/


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