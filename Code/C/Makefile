# Include for Matrix and GSL.
INC =   -I$(HOME)/RV-Project/Code/C_Examples/MyLib/Matrix/ \
	-I$(HOME)/RV-Project/Code/C_Examples/MyLib/RNG/ \
	-I$(HOME)/RV-Project/Code/C_Examples/MyLib/RandomVariates \
	-I$(HOME)/Code/include/

GLIB = -L$(HOME)/Code/lib

OPTS = -O2 -g -Wall 
## OPTS = -g -Wall -pedantic

## Note: pedantic kills eigen.

RINC = $(shell R CMD config --cppflags)
RLNK = $(shell R CMD config --ldflags)

test2: test2.cpp FFBS.cpp FFBS.h
	g++ -o test2 test2.cpp $(OPTS) $(INC) $(GLIB) -lgsl

test3: test3.cpp
	g++ -o test3 test3.cpp $(OPTS) $(INC) $(GLIB) -lgsl

test_pgpar : test_pgpar.cpp PolyaGamma.hpp PolyaGammaPar.hpp
	g++ -o test_pgpar test_pgpar.cpp $(OPTS) $(INC) $(GLIB) -fopenmp -lgsl -lblas -llapack 

test_logit : test_logit.cpp Logit.hpp
	g++ -o test_logit test_logit.cpp $(OPTS) $(INC) $(GLIB) -lgsl -lblas -llapack

GLogit.so : LogitWrapper.h LogitWrapper.cpp Logit.hpp
	g++ -o Logit.so LogitWrapper.cpp $(OPTS) $(INC) $(GLIB) -shared -fPIC -lgsl -lblas -llapack

RLogit.so : LogitWrapper.h LogitWrapper.cpp Logit.hpp
	g++ -o Logit.so LogitWrapper.cpp -DUSE_R $(OPTS) $(INC) $(RINC) $(RLNK) -shared -fPIC -lblas -llapack

hh_lambda : HHWrapper.hpp HHWrapper.cpp
	g++ -o hh_lambda.so HHWrapper.cpp -DUSE_R $(OPTS) $(INC) $(RINC) $(RLNK) -shared -fPIC -lblas -llapack

FSF_nmix.so : FSF_nmix.hpp FSF_nmix.cpp
	g++ -o FSF_nmix.so FSF_nmix.cpp -DUSE_R $(OPTS) $(INC) $(RINC) $(RLNK) -shared -fPIC -lblas -llapack

FFBS.so : FFBS.h FFBS.cpp
	g++ -o FFBS.so FFBS.cpp -DUSE_R $(OPTS) $(INC) $(RINC) $(RLNK) -shared -fPIC

################################################################################

hh_lambda.c : HHWrapper.hpp HHWrapper.cpp
	R CMD COMPILE CXXFLAGS="-DUSE_R $(OPTS) $(INC) $(RINC) -fPIC" HHWrapper.cpp

hh_lambda.o : HHWrapper.o
	R CMD SHLIB HHWrapper.o "-shared $(RLNK) -lblas -llapack -o hh_lambda.so" 

hh_lambda.c.2 : HHWrapper.hpp HHWrapper.cpp
	g++ -o HHWrapper.o -c HHWrapper.cpp -DUSE_R $(OPTS) $(INC) $(RINC) -fPIC

hh_lambda.o.2 : HHWrapper.o
	g++ -shared HHWrapper.o $(RLNK) -lblas -llapack -o HHWrapper.so

hh_lambda.gsl : HHWrapper.hpp HHWrapper.cpp
	g++ -o hh_lambda.so HHWrapper.cpp $(OPTS) $(INC) $(RINC) $(RLNK) $(GLIB) -shared -fPIC -lblas -llapack -lgsl