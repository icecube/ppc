ifdef I3_PORTS
dir	=	-I$(I3_PORTS)/include -L$(I3_PORTS)/lib
else
dir	=
endif
gsl	=	$(dir) -lgsl -lgslcblas

gcpp	=	$(CXX) -O2 #--fast-math

all:
	$(gcpp) bfr.cxx -o bfr $(gsl)
