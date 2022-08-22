ifdef I3_PORTS
dir	=	-I$(I3_PORTS)/include -L$(I3_PORTS)/lib
else
dir	=
endif
gsl	=	$(dir) -lgsl -lgslcblas

ifdef opcl
opts	=	-I$(opcl)/include -L$(opcl)/lib/x86_64 -lOpenCL
else
opts	=	-framework OpenCL
endif

gcpp	=	$(CXX) -O2 #--fast-math
nvcc	=	nvcc -O2 --use_fast_math

all:
	@echo "	make cpu:  compile against the  CPU   version in gpu directory"
	@echo "	make gpu:  compile against the  CUDA  version in gpu directory"
	@echo "	make ocl:  compile against the OpenCL version in ocl directory"

cpu:
	$(gcpp) llh.cxx -o llh $(gsl) ../gpu/ppc.o

gpu:
	$(nvcc) llh.cxx -o llh $(gsl) ../gpu/ppc.o

ocl:
	$(gcpp) llh.cxx -o llh $(gsl) ../ocl/ppc.o $(opts)

inv:
	c++ inv.cxx -o inv $(gsl)

ico:
	c++ ico.cxx -o ico
