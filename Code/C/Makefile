## PARAMETERS

# Include for Matrix and GSL.
UINC = -I./include -I./include/Matrix -I./include/RNG

GINC = -I$(HOME)/Code/include
GLIB = -L$(HOME)/Code/lib

RINC = $(shell R CMD config --cppflags)
RLNK = $(shell R CMD config --ldflags)

# USE_R = -DUSE_R
USE_R = 

# ## For BLAS and LAPACK -- I need to figure out how to find these libraries for any system.
# MKLLNK = -undefined dynamic_lookup -single_module -multiply_defined suppress \
#        -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRblas -lgfortran \
#        -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRlapack

VANLNK = -lblas -llapack

MATLNK = $(VANLNK)

ifndef USE_R
INC = $(UINC) $(GINC)
LNK = $(GLIB) $(MATLNK) -lgsl
else
INC = $(UINC) $(RINC)
LNK = $(MATLNK) $(RLNK) 
endif

OPTS = -O2 $(USE_R)
## OPTS = -g -Wall -pedantic $(USE_R)
## Note: pedantic kills eigen.

objects = LogitWrapper.so FSF_nmix.so FFBS.so hh_lambda.so

## TESTS -- set USE_R = ""

test2: test2.cpp FFBS.cpp FFBS.h
	g++ -o test2 test2.cpp $(OPTS) $(INC) $(LNK) 

test3: test3.cpp
	g++ -o test3 test3.cpp $(OPTS) $(INC) $(LNK)

test_pgpar : test_pgpar.cpp PolyaGamma.hpp PolyaGammaPar.hpp
	g++ -o test_pgpar test_pgpar.cpp $(OPTS) $(INC) -fopenmp $(LNK)

test_logit : test_logit.cpp Logit.hpp
	g++ -o test_logit test_logit.cpp $(OPTS) $(INC) $(LNK)

## OBJECTS

BayesLogit.so : $(objects)
	g++ -o BayesLogit.so $(objects)	$(OPTS) -shared -fPIC $(LNK)

LogitWrapper.so : LogitWrapper.hpp LogitWrapper.cpp Logit.hpp
	g++ -o LogitWrapper.so LogitWrapper.cpp $(OPTS) $(INC) -shared -fPIC $(LNK)

FSF_nmix.so : FSF_nmix.hpp FSF_nmix.cpp
	g++ -o FSF_nmix.so FSF_nmix.cpp $(OPTS) $(INC) -shared -fPIC $(LNK)

FFBS.so : FFBS.h FFBS.cpp
	g++ -o FFBS.so FFBS.cpp $(OPTS) $(INC) -shared -fPIC $(LNK)

hh_lambda.so : HHWrapper.hpp HHWrapper.cpp
	g++ -o hh_lambda.so HHWrapper.cpp $(OPTS) $(INC) -shared -fPIC $(LNK)

## CLEAN

clean: 
	rm *.so

cleandebug:
	rm -r *.dSYM

################################################################################
				 ## APPENDIX ##
################################################################################

## Link order matters.  First symbol found is used.