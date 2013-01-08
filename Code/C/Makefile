## PARAMETERS

# Include for Matrix and GSL.
UINC = -I./include -I./include/Matrix -I./include/RNG

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

OPTS = -O2 $(USE_R)
## OPTS = -g -Wall -pedantic $(USE_R)
## Note: pedantic kills eigen.

sobj = LogitWrapper.so FSF_nmix.so FFBS.so hh_lambda.so
obj = LogitWrapper.o FSF_nmix.o AR1.o CUBS.o CUBS_update.o FFBS.o hh_lambda.o RNG.o Matrix.o MatrixFrame.o $(DEP)

## TESTS -- set USE_R = ""

test2: test2.cpp FFBS.cpp FFBS.h
	g++ -o test2 test2.cpp $(OPTS) $(INC) $(LNK) 

test3: test3.cpp
	g++ -o test3 test3.cpp $(OPTS) $(INC) $(LNK)

test_pgpar : test_pgpar.cpp PolyaGamma.hpp PolyaGammaPar.hpp
	g++ -o test_pgpar BayesLogit.so test_pgpar.cpp $(OPTS) $(INC) -fopenmp $(LNK)

test_logit : test_logit.cpp Logit.hpp MatRNG.so
	g++ -o test_logit MatRNG.so test_logit.cpp $(OPTS) $(INC) $(LNK)

################################################################################

MatRNG.so : Matrix.o RNG.o
	g++ $(OPTS) $(INC) -o MatRNG.so Matrix.o MatrixFrame.o RNG.o $(DEP) -fPIC -shared $(LNK)

## Link gsl last so that we will use R libraries when compiling for R.
BayesLogit.so : $(obj)
	g++ -o BayesLogit.so $(obj) $(OPTS) -shared -fPIC $(LNK) $(GLIB) -lgsl

LogitWrapper.o : LogitWrapper.hpp LogitWrapper.cpp Logit.hpp
	g++ $(OPTS) $(INC) -o LogitWrapper.o -c LogitWrapper.cpp

FSF_nmix.o : FSF_nmix.hpp FSF_nmix.cpp
	g++ $(OPTS) $(INC) -o FSF_nmix.o -c FSF_nmix.cpp 

AR1.o : AR1.h AR1.cpp
	g++ $(OPTS) $(INC) -o AR1.o -c AR1.cpp

CUBS.o : CUBS.h CUBS.cpp
	g++ $(OPTS) $(INC) $(GINC) -o CUBS.o -c CUBS.cpp

CUBS_update.o : CUBS_update.cpp CUBS_update.h
	g++ $(OPTS) $(INC) $(GINC) -o CUBS_update.o -c CUBS_update.cpp

FFBS.o : FFBS.h FFBS.cpp
	g++ $(OPTS) $(INC) -o FFBS.o -c FFBS.cpp 

hh_lambda.o : HHWrapper.hpp HHWrapper.cpp
	g++ $(OPTS) $(INC) -o hh_lambda.o -c HHWrapper.cpp 

Matrix.o : MatrixFrame.o 
	g++ $(OPTS) $(INC) -c include/Matrix/Matrix.cpp -o Matrix.o

MatrixFrame.o :
	g++ $(OPTS) $(INC) -c include/Matrix/MatrixFrame.cpp -o MatrixFrame.o

RNG.o : $(DEP)
	g++ $(OPTS) $(INC) -c include/RNG/RNG.cpp -o RNG.o 

GRNG.o: 
	g++ $(INC) $(OPTs) -c include/RNG/GRNG.cpp -o GRNG.o

RRNG.o: 
	g++ $(INC) $(OPTs) -c include/RNG/RRNG.cpp -o RRNG.o

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

## You can created a shared object with any name and use it in R via dyn.load.
## However, R will created a shared object with the same name as the package.
## Thus you might as well do the same thing so that your wrapper functions use
## the correct parameters.