## PARAMETERS

# Include for Matrix and GSL.
UINC = -I./include -I./include/Matrix -I./include/RNG

# If you have manually installed GSL, set these to the appropriate directories.
GINC = -I$(HOME)/Code/include
GLIB = -L$(HOME)/Code/lib

RINC = $(shell R CMD config --cppflags)
RLNK = $(shell R CMD config --ldflags)

# make USE_R=-DUSE_R whatever
# USE_R = -DUSE_R
# USE_R = 

# ## For BLAS and LAPACK -- I need to figure out how to find these libraries for any system.
# MKLLNK = -undefined dynamic_lookup -single_module -multiply_defined suppress \
#        -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRblas -lgfortran \
#        -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRlapack

VANLNK = -lblas -llapack

MATLNK = $(VANLNK)

ifndef USE_R
INC = $(UINC) $(GINC)
LNK = $(GLIB) $(MATLNK) -lgsl
DEP = GRNG.o
else
INC = $(UINC) $(RINC)
LNK = $(MATLNK) $(RLNK) 
DEP = RRNG.o
endif

OPTS = -O2 $(USE_R) -g -DNDEBUG
## OPTS = -g -Wall -std=c++98  $(USE_R) ## -Wshadow
## OPTS = -g -Wall -pedantic $(USE_R)
## Note: pedantic kills eigen.

obj.neigen = LogitWrapper.o FSF_nmix.o PolyaGamma.o PolyaGammaAlt.o PolyaGammaSP.o InvertY.o InvertY2.o hh_lambda.o RNG.o Matrix.o MatrixFrame.o $(DEP)
obj.weigen = DynExpFamMH.o  FFBS.o AR1.o CUBS.o CUBS_update.o 

# Include Eigen stuff.
# obj = $(obj.neigen) $(obj.weigen)
# Don't include.
obj = $(obj.neigen)

sobj = LogitWrapper.so FSF_nmix.so FFBS.so hh_lambda.so

## TESTS -- set USE_R = ""

test2: test2.cpp FFBS.cpp FFBS.h
	g++ -o test2 test2.cpp $(OPTS) $(INC) $(LNK) 

test3: test3.cpp
	g++ -o test3 test3.cpp $(OPTS) $(INC) $(LNK)

test_hybrid_par : test_hybrid_par.cpp PolyaGamma.o PolyaGammaAlt.o PolyaGammaSP.o InvertY.o ParallelWrapper.o RNG.o
	g++ -o test_hybrid_par test_hybrid_par.cpp PolyaGamma.o PolyaGammaAlt.o PolyaGammaSP.o InvertY.o ParallelWrapper.o RNG.o GRNG.o $(OPTS) $(INC) -fopenmp $(LNK)

test_pgpar : test_pgpar.cpp PolyaGamma.h PolyaGammaPar.h
	g++ -o test_pgpar BayesLogit.so test_pgpar.cpp $(OPTS) $(INC) -fopenmp $(LNK)

test_logit : test_logit.cpp Logit.hpp MatRNG.so PolyaGamma.o PolyaGammaAlt.o PolyaGammaSP.o InvertY.o InvertY2.o
	g++ -o test_logit MatRNG.so PolyaGamma.o PolyaGammaAlt.o PolyaGammaSP.o InvertY.o InvertY2.o test_logit.cpp $(OPTS) $(INC) $(LNK)

test_cubs : CUBS_update.cpp CUBS_update.h 
	g++ -DCUBS_UPDATE_MAIN -o test_cubs CUBS_update.cpp $(OPTS) $(INC) $(LNK)

test_pgalt : test_pgalt.c PolyaGammaAlt.hpp BayesLogit.so
	g++ -o test_pgalt BayesLogit.so test_pgalt.c $(OPTS) $(INC) $(LNK)

test_pgsp : test_pgsp.c PolyaGammaSP.hpp BayesLogit.so 
	g++ -o test_pgsp BayesLogit.so test_pgsp.c $(OPTS) $(INC) $(LNK)

test_dynmh : test_dynmh.cpp DynExpFamMH.o MatRNG.so
	g++ -o test_dynmh MatRNG.so DynExpFamMH.o test_dynmh.cpp $(OPTS) $(INC) $(LNK)

################################################################################

## Link gsl last so that we will use R libraries when compiling for R.
BayesLogit.so : $(obj)
	g++ -o BayesLogit.so $(obj) $(OPTS) -shared -fPIC $(LNK) $(GLIB) -lgsl

MatRNG.so : Matrix.o RNG.o
	g++ $(OPTS) $(INC) -o MatRNG.so Matrix.o MatrixFrame.o RNG.o $(DEP) -fPIC -shared $(LNK)

ParallelWrapper.o : ParallelWrapper.cpp ParallelWrapper.h
	g++ $(OPTS) $(INC) -o ParallelWrapper.o -c ParallelWrapper.cpp -fopenmp

LogitWrapper.o : LogitWrapper.h LogitWrapper.cpp Logit.hpp PolyaGamma.o PolyaGammaAlt.o PolyaGammaSP.o
	g++ $(OPTS) $(INC) $(GINC) -fPIC -o LogitWrapper.o -c LogitWrapper.cpp

FSF_nmix.o : FSF_nmix.hpp FSF_nmix.cpp
	g++ $(OPTS) $(INC) -fPIC -o FSF_nmix.o -c FSF_nmix.cpp 

DynExpFamMH.o : DynExpFamMH.hpp DynExpFamMH.cpp
	g++ $(OPTS) $(INC) $(GINC) -fPIC -o DynExpFamMH.o -c DynExpFamMH.cpp

AR1.o : AR1.h AR1.cpp
	g++ $(OPTS) $(INC) $(GINC) -fPIC -o AR1.o -c AR1.cpp

CUBS.o : CUBS.h CUBS.cpp
	g++ $(OPTS) $(INC) $(GINC) -fPIC -o CUBS.o -c CUBS.cpp

CUBS_update.o : CUBS_update.cpp CUBS_update.h
	g++ $(OPTS) $(INC) $(GINC) -fPIC -o CUBS_update.o -c CUBS_update.cpp

PolyaGammaSP.o : PolyaGammaSP.h PolyaGammaSP.cpp
	g++ $(OPTS) $(INC) -fPIC -o PolyaGammaSP.o -c PolyaGammaSP.cpp

PolyaGammaAlt.o : PolyaGammaAlt.h PolyaGammaAlt.cpp
	g++ $(OPTS) $(INC) -fPIC -o PolyaGammaAlt.o -c PolyaGammaAlt.cpp

PolyaGamma.o : PolyaGamma.h PolyaGamma.cpp
	g++ $(OPTS) $(INC) -fPIC -o PolyaGamma.o -c PolyaGamma.cpp

InvertY2.o : InvertY2.hpp InvertY2.cpp 
	g++ $(OPTS) $(INC) $(GINC) -fPIC -o InvertY2.o -c InvertY2.cpp

InvertY.o : InvertY.hpp InvertY.cpp
	g++ $(OPTS) $(INC) -fPIC -o InvertY.o -c InvertY.cpp

FFBS.o : FFBS.h FFBS.cpp
	g++ $(OPTS) $(INC) $(GINC) -fPIC -o FFBS.o -c FFBS.cpp 

hh_lambda.o : HHWrapper.hpp HHWrapper.cpp
	g++ $(OPTS) $(INC) -fPIC -o hh_lambda.o -c HHWrapper.cpp 

Matrix.o : MatrixFrame.o 
	g++ $(OPTS) $(INC) -fPIC -c include/Matrix/Matrix.cpp -o Matrix.o

MatrixFrame.o :
	g++ $(OPTS) $(INC) -fPIC -c include/Matrix/MatrixFrame.cpp -o MatrixFrame.o

RNG.o : $(DEP)
	g++ $(OPTS) $(INC) -fPIC -c include/RNG/RNG.cpp -o RNG.o 

GRNG.o: 
	g++ $(INC) $(OPTS) -fPIC -c include/RNG/GRNG.cpp -o GRNG.o

RRNG.o: 
	g++ $(INC) $(OPTS) -fPIC -c include/RNG/RRNG.cpp -o RRNG.o

## CLEAN

clean: 
	rm *.o

cleandebug:
	rm -r *.dSYM

################################################################################
				 ## APPENDIX ##
################################################################################

## Link order matters.  First symbol found is used.

## I believe the way that I am making BayesLogit.so is now more in line with how
## R does things.

## You can create a shared object with any name and use it in R via dyn.load.
## However, R will create a shared object with the same name as the package.
## Thus you might as well do the same thing so that your wrapper functions use
## the correct parameters.